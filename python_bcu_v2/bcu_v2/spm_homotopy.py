"""SPM 联合平衡方程的连续同伦诊断器。

该模块用于独立检查 MGP 种子能否连续连接到物理 SPM 联合平衡根。它不读取
MATLAB 的 ``E_critical``，也不改变 ``spm_cuep`` 的生产接口；每个同伦步都以
上一步的完整网络状态作为初值，并把失败保留为结构化结果。
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping

import numpy as np
from scipy.optimize import least_squares

from . import spm_energy
from .spm_cuep import _project_coi
from .spm_joint_audit import _joint_residual


@dataclass
class SpmHomotopyResult:
    """同伦追踪的最终状态和逐步诊断。"""

    start_delta: np.ndarray
    delta_gen: np.ndarray
    theta_net: np.ndarray
    voltage_net: np.ndarray
    lambda_reached: float
    steps: list[dict[str, Any]]
    residual_norm: float
    network_residual: float
    equilibrium_residual: float
    branch_continuity_error: float
    e_critical: float
    equilibrium_type: str
    converged: bool
    used_external_ecritical: bool
    exit_reason: str

    def as_dict(self) -> dict[str, Any]:
        """返回可直接写入 JSON 的副本。"""

        payload = asdict(self)
        for key in ("start_delta", "delta_gen", "theta_net", "voltage_net"):
            payload[key] = np.asarray(payload[key], dtype=float).tolist()
        return payload


def _coerce_start(static, value: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    """严格读取完整物理网络种子，不填充缺失值。"""

    if not isinstance(value, Mapping):
        raise ValueError("start_state must be a mapping")
    required = ("delta_gen", "theta_net", "voltage_net")
    if any(name not in value for name in required):
        raise ValueError("start_state must contain delta_gen, theta_net and voltage_net")
    masses = np.asarray(static.preset.m, dtype=float).reshape(-1)
    ngen = int(masses.size)
    yfull = np.asarray(static.postfault.metadata.get("yfull_mod", static.postfault.yfull),
                       dtype=complex)
    nnet = int(yfull.shape[0] - ngen)
    delta = np.asarray(value["delta_gen"], dtype=float).reshape(-1)
    theta = np.asarray(value["theta_net"], dtype=float).reshape(-1)
    voltage = np.asarray(value["voltage_net"], dtype=float).reshape(-1)
    if delta.size == ngen - 1:
        delta = np.r_[-np.dot(masses[1:], delta) / masses[0], delta]
    if delta.size != ngen or theta.size != nnet or voltage.size != nnet:
        raise ValueError("start_state arrays have incompatible dimensions")
    if not all(np.all(np.isfinite(item)) for item in (delta, theta, voltage)):
        raise ValueError("start_state arrays must be finite")
    if not np.all(voltage > 1e-4):
        raise ValueError("start_state voltage_net must be positive")
    return _project_coi(delta, masses), np.r_[theta, voltage]


def _diagnose_state(static, delta: np.ndarray, network: np.ndarray,
                    sep: np.ndarray, sep_network: np.ndarray) -> dict[str, Any]:
    """计算单个同伦状态的残差、分类和自足势能。"""

    preset, post = static.preset, static.postfault
    yfull = np.asarray(post.metadata.get("yfull_mod", post.yfull), dtype=complex)
    ngen = int(np.asarray(preset.m).size)
    nnet = int(yfull.shape[0] - ngen)
    residual = _joint_residual(np.r_[delta, network], preset, yfull)
    network_residual = float(np.linalg.norm(
        spm_energy.spm_network_residual(network, delta, yfull, preset.epu)
    ))
    pe = spm_energy.spm_generator_power(
        delta, network[:nnet], network[nnet:], yfull, np.asarray(preset.epu, dtype=float)
    )
    mismatch = np.asarray(preset.pmpu, dtype=float) - pe
    damping = np.asarray(preset.d, dtype=float)
    omega = float(np.sum(mismatch) / np.sum(damping))
    equilibrium_residual = float(np.linalg.norm(mismatch - damping * omega))
    energy = spm_energy.spm_potential_energy(
        preset, post, yfull, sep, sep_network[:nnet], sep_network[nnet:],
        delta, network[:nnet], network[nnet:],
    )
    try:
        from .cuep import _reduced_jacobian_eig

        positive = int(np.sum(_reduced_jacobian_eig(delta, post.yred, preset) > 1e-6))
        equilibrium_type = "SEP" if positive == 0 else f"type-{positive}"
    except Exception:  # noqa: BLE001
        equilibrium_type = "UNVERIFIED"
    return {
        "residual_norm": float(np.linalg.norm(residual)),
        "network_residual": network_residual,
        "equilibrium_residual": equilibrium_residual,
        "e_critical": float(np.sum(energy)),
        "equilibrium_type": equilibrium_type,
        "voltage_min": float(np.min(network[nnet:])),
    }


def continue_spm_joint_homotopy(
    static,
    start_state: Mapping[str, Any],
    *,
    steps: int = 100,
    root_tol: float = 1e-10,
    residual_tol: float = 1e-8,
    max_nfev: int = 3000,
) -> SpmHomotopyResult:
    """从 MGP 状态连续追踪到 SPM 联合平衡根。

    令 ``F(x)`` 为发电机—网络联合残差，起点为 ``x0``，每个参数步求解
    ``F(x) - (1-lambda)F(x0) = 0``。在 ``lambda=0`` 时起点严格满足同伦方程，
    随后使用上一点的完整网络角/电压作为 warm-start。这个诊断只证明所登记
    起点上的连续分支，不宣称周期状态空间全根穷尽。
    """

    if int(steps) <= 0:
        raise ValueError("steps must be positive")
    if float(root_tol) <= 0 or float(residual_tol) <= 0:
        raise ValueError("root_tol and residual_tol must be positive")
    if int(max_nfev) <= 0:
        raise ValueError("max_nfev must be positive")

    preset, post = static.preset, static.postfault
    yfull = np.asarray(post.metadata.get("yfull_mod", post.yfull), dtype=complex)
    masses = np.asarray(preset.m, dtype=float).reshape(-1)
    ngen = int(masses.size)
    nnet = int(yfull.shape[0] - ngen)
    start_delta, start_network = _coerce_start(static, start_state)
    sep = _project_coi(np.asarray(post.sep_delta, dtype=float), masses)
    sep_network, sep_ok, sep_residual = spm_energy.solve_spm_network(
        sep, yfull, np.asarray(preset.epu, dtype=float), tol=1e-12,
    )
    if not sep_ok:
        raise RuntimeError(f"SEP network solve failed: residual={sep_residual:g}")

    current = np.r_[start_delta, start_network]
    initial_residual = _joint_residual(current, preset, yfull)
    records: list[dict[str, Any]] = []
    state_jumps: list[float] = []
    last_lambda = 0.0
    failure_reason = ""
    for index, lam in enumerate(np.linspace(1.0 / int(steps), 1.0, int(steps)), start=1):
        target = (1.0 - float(lam)) * initial_residual
        previous = current.copy()
        try:
            solution = least_squares(
                lambda value: _joint_residual(value, preset, yfull) - target,
                current, xtol=float(root_tol), ftol=float(root_tol), gtol=float(root_tol),
                max_nfev=int(max_nfev),
            )
            raw_delta = _project_coi(solution.x[:ngen], masses)
            network = np.asarray(solution.x[ngen:], dtype=float)
            if not (np.all(np.isfinite(raw_delta)) and np.all(np.isfinite(network))):
                raise ValueError("homotopy produced a non-finite state")
            if not np.all(network[nnet:] > 1e-4):
                raise ValueError("homotopy reached a nonphysical voltage branch")
            current = np.r_[raw_delta, network]
            diagnostics = _diagnose_state(static, raw_delta, network, sep, sep_network)
            jump = float(np.linalg.norm(current - previous))
            state_jumps.append(jump)
            diagnostics.update({"index": index, "lambda": float(lam),
                                "solver_success": bool(solution.success),
                                "state_jump": jump})
            records.append(diagnostics)
            last_lambda = float(lam)
            if not solution.success:
                failure_reason = f"homotopy solver stopped at lambda={lam:.6g}"
                break
        except Exception as exc:  # noqa: BLE001
            failure_reason = f"homotopy step {index} failed: {type(exc).__name__}: {exc}"
            break

    delta = _project_coi(current[:ngen], masses)
    network = np.asarray(current[ngen:], dtype=float)
    diagnostics = _diagnose_state(static, delta, network, sep, sep_network)
    reached = bool(last_lambda >= 1.0 - 1e-12)
    converged = bool(
        reached
        and diagnostics["network_residual"] < float(residual_tol)
        and diagnostics["equilibrium_residual"] < float(residual_tol)
        and diagnostics["equilibrium_type"].startswith("type-")
        and np.all(network[nnet:] > 1e-4)
    )
    if converged:
        reason = "joint homotopy reached physical type-1 root"
    elif failure_reason:
        reason = failure_reason
    else:
        reason = "joint homotopy finished without acceptance"
    return SpmHomotopyResult(
        start_delta=start_delta,
        delta_gen=delta,
        theta_net=network[:nnet],
        voltage_net=network[nnet:],
        lambda_reached=float(last_lambda),
        steps=records,
        residual_norm=float(diagnostics["residual_norm"]),
        network_residual=float(diagnostics["network_residual"]),
        equilibrium_residual=float(diagnostics["equilibrium_residual"]),
        branch_continuity_error=float(max(state_jumps, default=0.0)),
        e_critical=float(diagnostics["e_critical"]),
        equilibrium_type=str(diagnostics["equilibrium_type"]),
        converged=converged,
        used_external_ecritical=False,
        exit_reason=reason,
    )


__all__ = ["SpmHomotopyResult", "continue_spm_joint_homotopy"]
