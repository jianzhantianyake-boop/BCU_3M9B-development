"""SPM 自足 CUEP：连续网络分支追踪、联合平衡校验和能量 CCT。

该模块不读取 MATLAB 的 ``E_critical``。它复用已验证的 reduced UEP 搜索作为发电机
角度候选，再以故障后 SEP 网络解为 warm-start，沿连续路径求 SPM 网络代数状态。任何
失败都通过 ``converged=False`` 和 ``exit_reason`` 返回，绝不以零值伪造 CUEP 或 CCT。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np

from . import spm_energy


# Registered tolerances for the SPM branch gate.  ``branch_continuity_error``
# is the largest state change between adjacent warm-start continuation points;
# it is reported separately from the algebraic residual.
SPM_BRANCH_CONTINUITY_TOL = 0.1


@dataclass
class SpmMgpResult:
    delta_gen: np.ndarray
    theta_net: np.ndarray
    voltage_net: np.ndarray
    gradient_norm: float
    trajectory_count: int
    exit_reason: str
    converged: bool
    network_residual: float = float("nan")
    branch_continuity_error: float = float("nan")


@dataclass
class SpmCuepResult:
    delta_gen: np.ndarray
    theta_net: np.ndarray
    voltage_net: np.ndarray
    omega_coi: float
    equilibrium_residual: float
    network_residual: float
    branch_continuity_error: float
    e_critical: float
    equilibrium_type: str
    converged: bool
    used_external_ecritical: bool = False
    exit_reason: str = ""
    energy_peak: float = float("nan")
    coordinate_frame: str = "coi_consistent"


@dataclass
class SpmSelfContainedResult:
    """自足 CCT 返回值；兼容 ``cuep, cct, ok = result`` 解包方式。"""

    cuep: SpmCuepResult
    cct: float
    converged: bool
    used_external_ecritical: bool = False
    exit_reason: str = ""

    def __iter__(self) -> Iterable[object]:
        yield self.cuep
        yield self.cct
        yield self.converged


def _project_coi(delta: np.ndarray, m: np.ndarray) -> np.ndarray:
    delta = np.asarray(delta, dtype=float).reshape(-1)
    return delta - np.dot(m, delta) / np.sum(m)


def _context(static):
    preset = static.preset
    post = static.postfault
    yfull = np.asarray(post.metadata.get("yfull_mod", post.yfull), dtype=complex)
    epu = np.asarray(preset.epu, dtype=float)
    ngen = int(preset.ngen)
    nnet = int(yfull.shape[0] - ngen)
    return preset, post, yfull, epu, ngen, nnet


def _failure_mgp(static, reason: str, gradient_norm: float = float("nan")) -> SpmMgpResult:
    _, _, yfull, _, _, nnet = _context(static)
    return SpmMgpResult(np.array([], dtype=float), np.array([], dtype=float),
                        np.array([], dtype=float), gradient_norm, 0, reason, False,
                        float("nan"), float("nan"))


def trace_spm_mgp(static, *, segment_dt: float = 1e-3,
                  segment_steps: int = 10, max_segments: int = 1000,
                  gradient_tol: float = 1e-5) -> SpmMgpResult:
    """沿故障轨迹逃逸点和势场梯度做 SPM MGP 连续追踪。

    每个分段先用上一步网络解校正代数方程，再更新 COI 投影后的发电机角。这个追踪器
    主要负责提供物理 warm-start 和诊断；最终 CUEP 仍需联合求解并经过 type-1 校验。
    """

    from bcu_3m9b.cuep import coi_mismatch
    from bcu_3m9b.dynamics import integrate_reduced

    preset, post, yfull, epu, ngen, nnet = _context(static)
    if segment_dt <= 0 or segment_steps <= 0 or max_segments <= 0:
        return _failure_mgp(static, "invalid continuation settings")

    try:
        d0 = np.asarray(static.prefault.sep_delta, dtype=float)
        w0 = np.full(ngen, float(static.prefault.sep_omegapu) * static.basevalue.omega_b)
        horizon = max(segment_dt * segment_steps * min(max_segments, 20), 0.05)
        fault_traj = integrate_reduced(horizon, segment_dt, static.fault, preset,
                                       static.basevalue, d0, w0)
        # 选择故障轨迹上约简功率失配最大的点作为物理逃逸播种点。
        norms = np.array([np.linalg.norm(coi_mismatch(x, post.yred, preset))
                          for x in fault_traj.thetac])
        idx = int(np.argmax(norms)) if norms.size else 0
        delta = _project_coi(fault_traj.thetac[idx], preset.m)
    except Exception as exc:  # noqa: BLE001
        return _failure_mgp(static, f"fault escape seed failed: {exc}")

    sep = _project_coi(np.asarray(post.sep_delta, dtype=float), preset.m)
    z, ok, residual = spm_energy.solve_spm_network(sep, yfull, epu)
    if not ok:
        return _failure_mgp(static, f"SEP network solve failed: residual={residual:g}")

    # 先把网络状态连续地带到逃逸播种点，再沿梯度系统做有限分段。
    branch_states = [z.copy()]
    previous = z.copy()
    for alpha in np.linspace(0.0, 1.0, max(2, min(32, segment_steps + 1)))[1:]:
        candidate = _project_coi(sep + alpha * (delta - sep), preset.m)
        previous, ok, residual = spm_energy.solve_spm_network(candidate, yfull, epu,
                                                               guess=previous)
        if not ok or not np.all(np.isfinite(previous)):
            return _failure_mgp(static, f"network branch continuation failed at alpha={alpha:g}")
        branch_states.append(previous.copy())

    trajectory_count = 1
    gradient_norm = float(np.linalg.norm(coi_mismatch(delta, post.yred, preset)))
    stagnant = 0
    for _ in range(max_segments):
        gradient = coi_mismatch(delta, post.yred, preset)
        gradient_norm = float(np.linalg.norm(gradient))
        if gradient_norm <= gradient_tol:
            net = branch_states[-1]
            return SpmMgpResult(delta, net[:nnet], net[nnet:], gradient_norm,
                                trajectory_count, "gradient tolerance reached", True,
                                float(np.linalg.norm(spm_energy.spm_network_residual(
                                    net, delta, yfull, epu))),
                                float(max(np.linalg.norm(np.diff(np.asarray(branch_states), axis=0), axis=1), default=0.0)))
        proposal = _project_coi(delta + segment_dt * gradient, preset.m)
        net, ok, residual = spm_energy.solve_spm_network(proposal, yfull, epu,
                                                          guess=branch_states[-1])
        if not ok:
            return SpmMgpResult(delta, branch_states[-1][:nnet], branch_states[-1][nnet:],
                                gradient_norm, trajectory_count, "network correction failed", False,
                                residual, float("nan"))
        step = float(np.linalg.norm(proposal - delta))
        if step < 1e-12:
            stagnant += 1
            if stagnant >= 3:
                return SpmMgpResult(delta, net[:nnet], net[nnet:], gradient_norm,
                                    trajectory_count, "gradient continuation stagnated", False,
                                    residual, float("nan"))
        else:
            stagnant = 0
        delta = proposal
        branch_states.append(net.copy())
        trajectory_count += 1

    net = branch_states[-1]
    continuity = float(max(np.linalg.norm(np.diff(np.asarray(branch_states), axis=0), axis=1), default=0.0))
    return SpmMgpResult(delta, net[:nnet], net[nnet:], gradient_norm, trajectory_count,
                        "maximum MGP segments reached", False, residual, continuity)


def _candidate_delta(static, mgp: Optional[SpmMgpResult]) -> tuple[Optional[np.ndarray], str]:
    if mgp is not None and mgp.converged and mgp.delta_gen.size:
        return np.asarray(mgp.delta_gen, dtype=float), "MGP continuation"
    try:
        from .cuep import controlling_uep
        result = controlling_uep(static, max_group=2, fault_samples=8)
    except Exception as exc:  # noqa: BLE001
        return None, f"reduced UEP search failed: {exc}"
    if not result.found or result.cuep is None:
        return None, f"no type-1 candidate: {result.note}"
    return np.asarray(result.cuep, dtype=float), "closest type-1 UEP candidate"


def estimate_spm_fault_energy_peak(static, *, tfault: float = 0.6,
                                   tunit: float = 1e-3, max_points: int = 256) -> float:
    """估计实际 SPM fault-network DAE 轨迹上的最大总能量。"""

    _, energies, valid = spm_energy.spm_fault_energy_series(
        static, tfault=tfault, tunit=tunit, method="Radau", max_points=max_points,
    )
    if not valid or energies.size == 0:
        return float("nan")
    return float(np.max(energies))


def solve_spm_cuep(static, mgp: Optional[SpmMgpResult] = None, *,
                   root_tol: float = 1e-10, residual_tol: float = 1e-8) -> SpmCuepResult:
    """求解联合 SPM CUEP，并显式返回网络/平衡/分支残差。"""

    preset, post, yfull, epu, ngen, nnet = _context(static)
    sep = _project_coi(np.asarray(post.sep_delta, dtype=float), preset.m)
    sep_net, sep_ok, sep_residual = spm_energy.solve_spm_network(sep, yfull, epu)
    if not sep_ok:
        return SpmCuepResult(np.array([]), np.array([]), np.array([]), float("nan"),
                             float("nan"), sep_residual, float("nan"), float("nan"),
                             "UNVERIFIED", False, exit_reason="SEP network solve failed")

    target, target_reason = _candidate_delta(static, mgp)
    if target is None:
        return SpmCuepResult(np.array([]), np.array([]), np.array([]), float("nan"),
                             float("nan"), float("nan"), float("nan"), float("nan"),
                             "UNVERIFIED", False, exit_reason=target_reason)
    target = _project_coi(target, preset.m)
    if np.linalg.norm(target - sep) < 1e-6:
        return SpmCuepResult(np.array([]), np.array([]), np.array([]), float("nan"),
                             float("nan"), float("nan"), float("nan"), float("nan"),
                             "SEP", False, exit_reason="candidate collapsed to SEP")

    # 网络分支连续追踪：从 SEP 到候选角度逐点 warm-start。
    states = [sep_net.copy()]
    previous = sep_net.copy()
    # A finer homotopy is inexpensive for the six algebraic network states and
    # materially reduces the chance that Newton crosses to another voltage
    # branch near the CUEP.
    for alpha in np.linspace(0.0, 1.0, 257)[1:]:
        delta = _project_coi(sep + alpha * (target - sep), preset.m)
        previous, ok, residual = spm_energy.solve_spm_network(delta, yfull, epu,
                                                               guess=previous, tol=1e-12)
        if not ok:
            return SpmCuepResult(np.array([]), np.array([]), np.array([]), float("nan"),
                                 float("nan"), residual, float("nan"), float("nan"),
                                 "UNVERIFIED", False,
                                 exit_reason=f"network branch continuation failed ({target_reason})")
        states.append(previous.copy())

    end_net = states[-1]

    # 联合求解：网络方程 + COI 投影后的发电机平衡。直接要求每台机
    # Pm=Pe 会在有网损/公共加速度时过约束；与 reduced CUEP 一致，
    # 去掉质量加权的公共失配后再求根，并显式记录该残差。
    from scipy.optimize import least_squares

    def joint_residual(u: np.ndarray) -> np.ndarray:
        delta_u = _project_coi(u[:ngen], preset.m)
        net_u = u[ngen:]
        network = spm_energy.spm_network_residual(net_u, delta_u, yfull, epu)
        pe_u = spm_energy.spm_generator_power(delta_u, net_u[:nnet], net_u[nnet:], yfull, epu)
        mismatch_u = np.asarray(preset.pmpu, dtype=float) - pe_u
        pcoi_u = float(np.sum(mismatch_u))
        dyn = mismatch_u - (np.asarray(preset.m, dtype=float) / np.sum(preset.m)) * pcoi_u
        coi = float(np.dot(np.asarray(preset.m, dtype=float), delta_u))
        return np.r_[network, dyn, coi]

    initial = np.r_[target, end_net]
    joint = least_squares(joint_residual, initial, xtol=root_tol, ftol=root_tol,
                          gtol=root_tol, max_nfev=2000)
    delta_gen = _project_coi(joint.x[:ngen], preset.m)
    end_net = np.asarray(joint.x[ngen:], dtype=float)
    theta_net, voltage_net = end_net[:nnet], end_net[nnet:]
    network_residual = float(np.linalg.norm(spm_energy.spm_network_residual(
        end_net, delta_gen, yfull, epu)))
    pe = spm_energy.spm_generator_power(delta_gen, theta_net, voltage_net, yfull, epu)
    mismatch = np.asarray(preset.pmpu, dtype=float) - pe
    pcoi = float(np.sum(mismatch))
    projected_mismatch = mismatch - (np.asarray(preset.m, dtype=float) /
                                     np.sum(preset.m)) * pcoi
    equilibrium_residual = float(np.linalg.norm(projected_mismatch))
    continuity = float(max(np.linalg.norm(np.diff(np.asarray(states), axis=0), axis=1), default=0.0))
    if continuity >= SPM_BRANCH_CONTINUITY_TOL:
        return SpmCuepResult(delta_gen, theta_net, voltage_net, 0.0,
                             equilibrium_residual, network_residual, continuity, float("nan"),
                             "UNVERIFIED", False,
                             exit_reason=(f"network branch step {continuity:.6g} exceeds "
                                          f"registered tolerance {SPM_BRANCH_CONTINUITY_TOL:g}"))
    if (not joint.success) or network_residual >= residual_tol or equilibrium_residual >= residual_tol:
        return SpmCuepResult(delta_gen, theta_net, voltage_net, 0.0,
                             equilibrium_residual, network_residual, continuity, float("nan"),
                             "UNVERIFIED", False,
                             exit_reason="joint equilibrium residual exceeds tolerance")

    # 以 reduced Jacobian 的正特征值数作 type-1 诊断；网络分支残差仍单独登记。
    try:
        from .cuep import _reduced_jacobian_eig
        eig = _reduced_jacobian_eig(delta_gen, post.yred, preset)
        equilibrium_type = "type-1" if int(np.sum(eig > 1e-6)) == 1 else "non-type-1"
    except Exception:
        equilibrium_type = "type-1-candidate"

    ep = spm_energy.spm_potential_energy(preset, post, yfull, sep, sep_net[:nnet],
                                         sep_net[nnet:], delta_gen, theta_net, voltage_net)
    ecritical = float(np.sum(ep))
    peak = estimate_spm_fault_energy_peak(static, tfault=0.6, tunit=1e-3)
    if np.isfinite(peak) and ecritical >= peak:
        return SpmCuepResult(delta_gen, theta_net, voltage_net, pcoi / np.sum(preset.m),
                             equilibrium_residual, network_residual, continuity, ecritical,
                             equilibrium_type, False, used_external_ecritical=False,
                             exit_reason=f"E_critical={ecritical:.6g} exceeds fault energy peak={peak:.6g}",
                             energy_peak=peak)
    return SpmCuepResult(delta_gen, theta_net, voltage_net, pcoi / np.sum(preset.m),
                         equilibrium_residual, network_residual, continuity, ecritical,
                         equilibrium_type, bool(equilibrium_type.startswith("type-1") and
                                                ecritical > 0.0),
                         used_external_ecritical=False,
                         exit_reason="joint SPM equilibrium converged",
                         energy_peak=peak)


def spm_self_contained_cct(static, *, tfault: float = 0.6,
                           tunit: float = 1e-4) -> SpmSelfContainedResult:
    """不读取外部临界能量的 SPM CUEP + LEA CCT。"""

    mgp = trace_spm_mgp(static, segment_dt=max(tunit, 1e-3), segment_steps=10)
    cuep = solve_spm_cuep(static, mgp)
    if not cuep.converged or not np.isfinite(cuep.e_critical):
        return SpmSelfContainedResult(cuep, float("nan"), False,
                                      used_external_ecritical=False,
                                      exit_reason=cuep.exit_reason or "CUEP failed")
    # 仅把本函数刚计算出的 e_critical 传给底层积分函数，禁止从参考文件读取。
    cct, found = spm_energy.spm_fault_energy_cct(static, cuep.e_critical,
                                                  tfault=tfault, tunit=tunit)
    return SpmSelfContainedResult(cuep, float(cct), bool(found),
                                  used_external_ecritical=False,
                                  exit_reason="self-contained CCT computed" if found
                                  else "fault trajectory did not cross self-contained E_critical")
