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


def _solve_network_continuous(delta_target: np.ndarray, yfull: np.ndarray,
                              epu: np.ndarray, guess: np.ndarray,
                              delta_anchor: np.ndarray | None = None,
                              *, max_angle_step: float = 2e-4,
                              tol: float = 1e-11) -> tuple[np.ndarray, bool, float]:
    """Follow one algebraic branch with small angle homotopy steps.

    A direct ``hybr`` call can converge to the zero-voltage mathematical root
    when a ray approaches a fold.  MATLAB's AE Newton receives every 1e-3 ray
    point sequentially; splitting each Python step into smaller continuation
    points reproduces that warm-start invariant while retaining a structured
    failure result.
    """

    target = np.asarray(delta_target, dtype=float).reshape(-1)
    anchor = target if delta_anchor is None else np.asarray(delta_anchor, dtype=float).reshape(-1)
    distance = float(np.linalg.norm(target - anchor))
    count = max(1, int(np.ceil(distance / max(float(max_angle_step), 1e-8))))
    initial = np.asarray(guess, dtype=float).reshape(-1).copy()

    def run(step: float):
        n = max(1, int(np.ceil(distance / max(float(step), 1e-8))))
        current = initial.copy()
        residual = float("inf")
        for k in range(1, n + 1):
            point = anchor + (target - anchor) * (k / n)
            current, ok, residual = spm_energy.solve_spm_network_newton(
                point, yfull, epu, guess=current, tol=tol,
            )
            if not ok:
                # Keep a conservative fallback for cases where the analytic
                # Newton Jacobian is singular, but never accept a nonphysical
                # voltage root as success.
                current, ok, residual = spm_energy.solve_spm_network(
                    point, yfull, epu, guess=current, tol=tol,
                )
                if not ok:
                    return current, False, residual
        return current, True, residual

    current, ok, residual = run(max_angle_step)
    if not ok and max_angle_step > 1.1e-4:
        # Near a voltage fold the 2e-4 default can land on the zero-voltage
        # root.  Retry the whole short continuation from the original warm
        # start with half-size steps before declaring the branch unavailable.
        current, ok, residual = run(max_angle_step / 2.0)
    return current, ok, residual


def _spm_escape_seed(static, *, tfault: float = 0.5,
                     tunit: float = 1e-4) -> tuple[np.ndarray, np.ndarray, float]:
    """Find the SPM fault dot-product crossing used by MATLAB MGP."""

    from .spm_dae import remap_algebraic_state, simulate_spm_dae

    preset, post, yfull, epu, ngen, nnet = _context(static)
    d0 = np.asarray(static.prefault.sep_delta, dtype=float)
    omega0 = np.full(ngen, float(static.prefault.sep_omegapu) * static.basevalue.omega_b)
    pref_y = np.asarray(static.prefault.metadata.get("yfull_mod", static.prefault.yfull), dtype=complex)
    pref_z, pref_ok, pref_residual = spm_energy.solve_spm_network(d0, pref_y, epu)
    if not pref_ok:
        raise RuntimeError(f"prefault network residual={pref_residual:g}")
    fault_guess = remap_algebraic_state(pref_z, static.prefault, static.fault, ngen)
    traj = simulate_spm_dae(tfault, tunit, static.fault, preset, static.basevalue,
                            d0, omega0, method="Radau", algebraic_guess=fault_guess)
    if not traj.get("success", False) or traj["time"].size < 3:
        raise RuntimeError("strict fault DAE did not converge")
    sep = _project_coi(np.asarray(post.sep_delta, dtype=float), preset.m)
    sep_z, sep_ok, sep_residual = spm_energy.solve_spm_network(sep, yfull, epu)
    if not sep_ok:
        raise RuntimeError(f"postfault SEP residual={sep_residual:g}")
    previous = sep_z.copy()
    dot = np.full(traj["time"].size, np.nan, dtype=float)
    crossing_index: int | None = None
    for k, (dg, omega) in enumerate(zip(traj["delta"], traj["omega"])):
        previous, ok, residual = spm_energy.solve_spm_network(
            np.asarray(dg, dtype=float), yfull, epu, guess=previous, tol=1e-11,
        )
        if not ok:
            # After the first PEBS crossing the fault trajectory may approach
            # a low-voltage algebraic branch.  The exit seed is already fixed
            # by then, so do not let a later nonphysical point erase valid
            # evidence.
            if crossing_index is not None:
                break
            raise RuntimeError(f"postfault exit-point network residual={residual:g}")
        pe = spm_energy.spm_generator_power(
            dg, previous[:nnet], previous[nnet:], yfull, epu,
        )
        omega_c = np.asarray(omega, dtype=float) - np.dot(preset.m, omega) / np.sum(preset.m)
        dot[k] = float(np.dot(np.asarray(preset.pmpu) - pe, omega_c))
        if k >= 1 and dot[k - 1] < -1e-9 and dot[k] > 1e-9:
            crossing_index = k - 1
            break
    if crossing_index is None:
        raise RuntimeError("SPM fault dot-product crossing not found")
    idx = int(crossing_index)
    dg = _project_coi(np.asarray(traj["delta"][idx], dtype=float), preset.m)
    z, ok, residual = spm_energy.solve_spm_network(dg, yfull, epu, guess=previous, tol=1e-11)
    if not ok:
        raise RuntimeError(f"escape postfault network residual={residual:g}")
    return dg, z, float(traj["time"][idx])


def _spm_gradient(delta: np.ndarray, z: np.ndarray, preset,
                  yfull: np.ndarray) -> np.ndarray:
    ngen = int(np.asarray(delta).size)
    nnet = int(yfull.shape[0]) - ngen
    pe = spm_energy.spm_generator_power(
        delta, np.asarray(z)[:nnet], np.asarray(z)[nnet:], yfull,
        np.asarray(preset.epu, dtype=float),
    )
    mismatch = np.asarray(preset.pmpu, dtype=float) - pe
    m = np.asarray(preset.m, dtype=float)
    return mismatch - m * np.sum(mismatch) / np.sum(m)


def _spm_mgp_trajectory(start_delta: np.ndarray, start_z: np.ndarray,
                         static, *, dt: float, steps: int,
                         norm_tol: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, float]:
    """One MATLAB-compatible ten-point MGP trajectory (RK4 on the DAE manifold)."""

    preset, post, yfull, epu, ngen, nnet = _context(static)
    m = np.asarray(preset.m, dtype=float)
    d = np.asarray(start_delta, dtype=float).copy()
    z = np.asarray(start_z, dtype=float).copy()
    x = d[1:].copy()
    ds: list[np.ndarray] = []
    zs: list[np.ndarray] = []
    norms: list[float] = []
    previous_d = d.copy()

    def solve_at(point: np.ndarray, guess: np.ndarray, anchor: np.ndarray) -> np.ndarray:
        out, ok, residual = _solve_network_continuous(
            point, yfull, epu, guess, anchor, max_angle_step=1e-3,
        )
        if not ok:
            raise RuntimeError(f"network residual={residual:g}")
        return out

    for output in range(max(1, int(steps))):
        dg = _project_coi(np.r_[-np.dot(m[1:], x) / m[0], x], m)
        z = solve_at(dg, z, previous_d)
        ds.append(dg.copy())
        zs.append(z.copy())
        norms.append(float(np.linalg.norm(_spm_gradient(dg, z, preset, yfull))))
        previous_d = dg.copy()
        if output == steps - 1:
            break

        def rhs(xx: np.ndarray, guess: np.ndarray, anchor: np.ndarray):
            dd = _project_coi(np.r_[-np.dot(m[1:], xx) / m[0], xx], m)
            zz = solve_at(dd, guess, anchor)
            return _spm_gradient(dd, zz, preset, yfull)[1:] / m[1:], zz

        k1, z1 = rhs(x, z, previous_d)
        k2, z2 = rhs(x + 0.5 * dt * k1, z1, previous_d)
        k3, z3 = rhs(x + 0.5 * dt * k2, z2, previous_d)
        k4, z4 = rhs(x + dt * k3, z3, previous_d)
        x = x + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        z = z4
    norms_a = np.asarray(norms, dtype=float)
    found = -1
    norm_min = norms_a[0]
    for k in range(1, norms_a.size):
        if norms_a[k] < norm_min:
            norm_min = norms_a[k]
        if (k >= 2 and norms_a[k] - norms_a[k - 1] > norm_tol
                and abs(norms_a[k - 1] - norm_min) <= 1e-10
                and norm_min < 1e-1):
            found = k - 1
            break
    return np.asarray(ds), np.asarray(zs), norms_a, found, float(norm_min)


def _spm_ray_update(last_delta: np.ndarray, last_z: np.ndarray, static,
                    *, ray_step: float = 1e-3,
                    path_energy_cal: int = 20) -> tuple[np.ndarray, np.ndarray, bool, float]:
    """Search a MATLAB-style local maximum on the SEP-to-last-point ray."""

    preset, post, yfull, epu, ngen, nnet = _context(static)
    m = np.asarray(preset.m, dtype=float)
    sep = _project_coi(np.asarray(post.sep_delta, dtype=float), m)
    sep_z, ok, residual = spm_energy.solve_spm_network(sep, yfull, epu)
    if not ok:
        raise RuntimeError(f"SEP network residual={residual:g}")
    direction = np.asarray(last_delta, dtype=float) - sep
    distance = float(np.linalg.norm(direction))
    if distance < 1e-12:
        return np.asarray(last_delta), np.asarray(last_z), False, float("nan")
    direction /= distance
    max_steps = max(1, int(np.floor(2.0 * distance / ray_step)))
    prev_d, prev_z = sep.copy(), sep_z.copy()
    prev_e = float(np.sum(spm_energy.spm_potential_energy(
        preset, post, yfull, sep, sep_z[:nnet], sep_z[nnet:],
        sep, sep_z[:nnet], sep_z[nnet:], path_energy_cal=path_energy_cal)))
    cur_d = sep + ray_step * direction
    cur_z, ok, residual = _solve_network_continuous(
        cur_d, yfull, epu, prev_z, prev_d, max_angle_step=1e-3,
    )
    if not ok:
        return np.asarray(last_delta), np.asarray(last_z), False, float("nan")
    cur_e = float(np.sum(spm_energy.spm_potential_energy(
        preset, post, yfull, sep, sep_z[:nnet], sep_z[nnet:],
        cur_d, cur_z[:nnet], cur_z[nnet:], path_energy_cal=path_energy_cal)))
    values = [prev_e, cur_e]
    for _ in range(2, max_steps + 1):
        nxt_d = cur_d + ray_step * direction
        nxt_z, ok, residual = _solve_network_continuous(
            nxt_d, yfull, epu, cur_z, cur_d, max_angle_step=1e-3,
        )
        if not ok:
            return np.asarray(last_delta), np.asarray(last_z), False, float("nan")
        nxt_e = float(np.sum(spm_energy.spm_potential_energy(
            preset, post, yfull, sep, sep_z[:nnet], sep_z[nnet:],
            nxt_d, nxt_z[:nnet], nxt_z[nnet:], path_energy_cal=path_energy_cal)))
        values.append(nxt_e)
        if values[-2] > values[-3] and values[-2] > values[-1]:
            return cur_d, cur_z, True, values[-2]
        prev_d, prev_z, prev_e = cur_d, cur_z, cur_e
        cur_d, cur_z, cur_e = nxt_d, nxt_z, nxt_e
    return np.asarray(last_delta), np.asarray(last_z), False, float(values[-1])


def trace_spm_mgp(static, *, segment_dt: float = 1e-3,
                  segment_steps: int = 10, max_segments: int = 1000,
                  gradient_tol: float = 1e-5,
                  path_energy_cal: int = 20) -> SpmMgpResult:
    """按 MATLAB 外层逻辑追踪 SPM MGP，并返回结构化失败原因。"""

    preset, post, yfull, epu, ngen, nnet = _context(static)
    if segment_dt <= 0 or segment_steps <= 0 or max_segments <= 0:
        return _failure_mgp(static, "invalid continuation settings")
    try:
        delta, z, escape_time = _spm_escape_seed(static, tfault=0.5, tunit=1e-4)
    except Exception as exc:  # noqa: BLE001
        return _failure_mgp(static, f"fault escape seed failed: {exc}")
    sep = _project_coi(np.asarray(post.sep_delta, dtype=float), preset.m)
    branch_states: list[np.ndarray] = [z.copy()]
    trajectory_count = 0
    gradient_norm = float("nan")
    last_residual = float("nan")
    try:
        for _ in range(int(max_segments)):
            ds, zs, norms, found, norm_min = _spm_mgp_trajectory(
                delta, z, static, dt=segment_dt, steps=segment_steps,
                norm_tol=gradient_tol,
            )
            trajectory_count += 1
            branch_states.extend(list(zs))
            gradient_norm = float(norms[-1])
            if found >= 0:
                mgp_delta, mgp_z = ds[found], zs[found]
                residual = float(np.linalg.norm(spm_energy.spm_network_residual(
                    mgp_z, mgp_delta, yfull, epu)))
                continuity = float(max(np.linalg.norm(np.diff(np.asarray(branch_states), axis=0), axis=1), default=0.0))
                return SpmMgpResult(
                    mgp_delta, mgp_z[:nnet], mgp_z[nnet:], float(norms[found]),
                    trajectory_count, "MGP local minimum reached", True,
                    residual, continuity,
                )
            last_delta, last_z = ds[-1], zs[-1]
            delta_new, z_new, found_ray, _ = _spm_ray_update(
                last_delta, last_z, static, path_energy_cal=path_energy_cal,
            )
            # MATLAB only applies the repeated-status termination when the
            # ray updater actually found a local maximum.  A ray with no
            # local maximum simply seeds the next trajectory at its last
            # point; conflating the two would stop after trajectory 1.
            if not found_ray:
                delta, z = last_delta, last_z
                continue
            delta, z = delta_new, z_new
            shift = float(np.linalg.norm(delta - last_delta))
            if shift < 1e-3:
                residual = float(np.linalg.norm(spm_energy.spm_network_residual(
                    z, delta, yfull, epu)))
                continuity = float(max(np.linalg.norm(np.diff(np.asarray(branch_states), axis=0), axis=1), default=0.0))
                return SpmMgpResult(
                    delta, z[:nnet], z[nnet:], gradient_norm, trajectory_count,
                    "MGP ray update repeated", True, residual, continuity,
                )
            if shift > 0.5 * float(np.linalg.norm(last_delta - sep)):
                delta, z = last_delta, last_z
            last_residual = float(np.linalg.norm(spm_energy.spm_network_residual(
                z, delta, yfull, epu)))
    except Exception as exc:  # noqa: BLE001
        return SpmMgpResult(
            np.asarray(delta), np.asarray(z)[:nnet], np.asarray(z)[nnet:],
            gradient_norm, trajectory_count, f"MGP continuation failed: {exc}",
            False, last_residual, float("nan"),
        )
    continuity = float(max(np.linalg.norm(np.diff(np.asarray(branch_states), axis=0), axis=1), default=0.0))
    return SpmMgpResult(delta, z[:nnet], z[nnet:], gradient_norm,
                        trajectory_count, "maximum MGP trajectories reached", False,
                        last_residual, continuity)


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
        # At a common-speed equilibrium the residual power is absorbed by
        # damping: PCOI = sum(d) * omega_coi.  Expressing the generator
        # balance directly as mismatch - d*omega avoids the old, incorrect
        # inertia normalization (which changes the returned speed offset by
        # a factor of sum(d)/sum(m)).
        damping = np.asarray(preset.d, dtype=float)
        pcoi_u = float(np.sum(mismatch_u))
        omega_u = pcoi_u / float(np.sum(damping))
        dyn = mismatch_u - damping * omega_u
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
    damping = np.asarray(preset.d, dtype=float)
    omega_coi = pcoi / float(np.sum(damping))
    equilibrium_residual = float(np.linalg.norm(mismatch - damping * omega_coi))
    continuity = float(max(np.linalg.norm(np.diff(np.asarray(states), axis=0), axis=1), default=0.0))
    if continuity >= SPM_BRANCH_CONTINUITY_TOL:
        return SpmCuepResult(delta_gen, theta_net, voltage_net, omega_coi,
                             equilibrium_residual, network_residual, continuity, float("nan"),
                             "UNVERIFIED", False,
                             exit_reason=(f"network branch step {continuity:.6g} exceeds "
                                          f"registered tolerance {SPM_BRANCH_CONTINUITY_TOL:g}"))
    if (not joint.success) or network_residual >= residual_tol or equilibrium_residual >= residual_tol:
        return SpmCuepResult(delta_gen, theta_net, voltage_net, omega_coi,
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
    # MATLAB Fun_Cal_CCT_Energy_SPM scans the registered 0.5 s fault
    # trajectory.  Use that same window for the physical E_critical gate;
    # the public self-contained CCT interface may still integrate farther
    # when its caller requests a different ``tfault``.
    peak = estimate_spm_fault_energy_peak(static, tfault=0.5, tunit=1e-3)
    if np.isfinite(peak) and ecritical >= peak:
        return SpmCuepResult(delta_gen, theta_net, voltage_net, omega_coi,
                             equilibrium_residual, network_residual, continuity, ecritical,
                             equilibrium_type, False, used_external_ecritical=False,
                             exit_reason=f"E_critical={ecritical:.6g} exceeds fault energy peak={peak:.6g}",
                             energy_peak=peak)
    return SpmCuepResult(delta_gen, theta_net, voltage_net, omega_coi,
                         equilibrium_residual, network_residual, continuity, ecritical,
                         equilibrium_type, bool(equilibrium_type.startswith("type-1") and
                                                ecritical > 0.0),
                         used_external_ecritical=False,
                         exit_reason="joint SPM equilibrium converged",
                         energy_peak=peak)


def spm_self_contained_cct(static, *, tfault: float = 0.6,
                           tunit: float = 1e-4,
                           max_segments: int = 1000) -> SpmSelfContainedResult:
    """不读取外部临界能量的 SPM CUEP + LEA CCT。"""

    if max_segments <= 0:
        raise ValueError("max_segments must be positive")
    mgp = trace_spm_mgp(static, segment_dt=max(tunit, 1e-3), segment_steps=10,
                        max_segments=max_segments)
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
