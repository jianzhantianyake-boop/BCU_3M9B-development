"""可重复的 SPM 发电机—网络联合平衡根审计。

这个模块是诊断工具，不会把任何 MATLAB ``E_critical`` 当作输入。它使用与
``spm_cuep.solve_spm_cuep`` 相同的联合残差，记录有限初值搜索中每个启动点的
收敛、网络残差、平衡残差、正电压条件、type-1 分类、角度展开和势能。有限搜索
只能扩大证据范围，不能宣称穷尽周期状态空间的全部数学根。
"""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
from scipy.optimize import least_squares

from . import spm_energy
from .spm_cuep import _project_coi


def _joint_residual(u: np.ndarray, preset, yfull: np.ndarray) -> np.ndarray:
    """与 ``solve_spm_cuep`` 一致的发电机/网络联合残差。"""

    ngen = int(np.asarray(preset.m).size)
    delta = _project_coi(np.asarray(u[:ngen], dtype=float), np.asarray(preset.m, dtype=float))
    net = np.asarray(u[ngen:], dtype=float)
    nnet = int(yfull.shape[0] - ngen)
    network = spm_energy.spm_network_residual(net, delta, yfull, preset.epu)
    pe = spm_energy.spm_generator_power(
        delta, net[:nnet], net[nnet:], yfull, np.asarray(preset.epu, dtype=float)
    )
    mismatch = np.asarray(preset.pmpu, dtype=float) - pe
    damping = np.asarray(preset.d, dtype=float)
    omega = float(np.sum(mismatch) / np.sum(damping))
    dyn = mismatch - damping * omega
    coi = float(np.dot(np.asarray(preset.m, dtype=float), delta))
    return np.r_[network, dyn, coi]


def _unwrap_to_reference(delta: np.ndarray, reference: np.ndarray,
                         masses: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """选择与参考状态最近的 ``2π`` 角度展开，并重新投影到 COI。"""

    raw = np.asarray(delta, dtype=float).reshape(-1)
    ref = np.asarray(reference, dtype=float).reshape(raw.size)
    wraps = np.rint((ref - raw) / (2.0 * np.pi)).astype(int)
    unwrapped = _project_coi(raw + 2.0 * np.pi * wraps, masses)
    return unwrapped, [int(item) for item in wraps]


def _network_branch_to(static, delta: np.ndarray, sep_delta: np.ndarray,
                       sep_state: np.ndarray, *, segments: int = 64) -> tuple[np.ndarray, float, bool]:
    """从 SEP 网络状态连续追踪到目标角，返回末态、最大步长和成功标志。"""

    preset, post = static.preset, static.postfault
    yfull = np.asarray(post.metadata.get("yfull_mod", post.yfull), dtype=complex)
    current = np.asarray(sep_state, dtype=float).copy()
    max_step = 0.0
    for alpha in np.linspace(1.0 / segments, 1.0, segments):
        point = (1.0 - alpha) * np.asarray(sep_delta, dtype=float) + alpha * np.asarray(delta, dtype=float)
        state, ok, _ = spm_energy.solve_spm_network(
            point, yfull, np.asarray(preset.epu, dtype=float), guess=current, tol=1e-11,
        )
        if not ok or not np.all(np.isfinite(state)):
            return current, float("inf"), False
        max_step = max(max_step, float(np.linalg.norm(state - current)))
        current = np.asarray(state, dtype=float)
    return current, max_step, True


def _classify(delta: np.ndarray, static) -> str:
    try:
        from .cuep import _reduced_jacobian_eig

        eig = _reduced_jacobian_eig(delta, static.postfault.yred, static.preset)
        npos = int(np.sum(np.asarray(eig) > 1e-6))
        return "SEP" if npos == 0 else ("type-1" if npos == 1 else f"type-{npos}")
    except Exception:  # noqa: BLE001
        return "unknown"


def _default_seeds(sep: np.ndarray, masses: np.ndarray, count: int,
                   random_seed: int) -> list[np.ndarray]:
    rng = np.random.default_rng(int(random_seed))
    seeds = [np.asarray(sep, dtype=float).copy()]
    for _ in range(max(0, int(count) - 1)):
        # The scale is deliberately broad enough to probe wrapped and distinct
        # generator-angle basins, while the deterministic seed makes the audit
        # reproducible.  Network states are always warm-started from SEP.
        candidate = np.asarray(sep, dtype=float) + rng.normal(0.0, 2.0, size=sep.size)
        seeds.append(_project_coi(candidate, masses))
    return seeds


def _coerce_seed_state(value, ngen: int, nnet: int,
                       masses: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Validate and normalize one explicit ``(delta, theta, voltage)`` seed.

    JSON reports use a mapping with three named arrays.  A flat numeric vector
    of length ``ngen + 2*nnet`` is accepted as a small convenience for callers
    constructing seeds programmatically.  No missing network values are
    filled in: an explicit state is either complete and finite or rejected.
    """

    if isinstance(value, Mapping):
        required = ("delta_gen", "theta_net", "voltage_net")
        if any(key not in value for key in required):
            raise ValueError("seed_states mapping must contain delta_gen, theta_net and voltage_net")
        delta = np.asarray(value["delta_gen"], dtype=float).reshape(-1)
        theta = np.asarray(value["theta_net"], dtype=float).reshape(-1)
        voltage = np.asarray(value["voltage_net"], dtype=float).reshape(-1)
        if theta.size != nnet or voltage.size != nnet:
            raise ValueError("seed_states network arrays have an incompatible length")
        state = np.r_[theta, voltage]
    else:
        flat = np.asarray(value, dtype=float).reshape(-1)
        expected = ngen + 2 * nnet
        if flat.size != expected:
            raise ValueError(f"flat seed state must contain {expected} values")
        delta, state = flat[:ngen], flat[ngen:]
    if delta.size == ngen - 1:
        delta = np.r_[-np.dot(masses[1:], delta) / masses[0], delta]
    if delta.size != ngen or not np.all(np.isfinite(delta)):
        raise ValueError("seed_states contains an incompatible or non-finite angle vector")
    if not np.all(np.isfinite(state)):
        raise ValueError("seed_states contains non-finite network values")
    return _project_coi(delta, masses), np.asarray(state, dtype=float)


def audit_spm_joint_roots(static, *, max_starts: int = 152,
                          random_seed: int = 20260902,
                          seed_deltas: Sequence[Sequence[float]] | None = None,
                          seed_states: Sequence[Mapping | Sequence[float]] | None = None,
                          root_tol: float = 1e-10,
                          residual_tol: float = 1e-8) -> dict:
    """对有限初值集合执行 SPM 联合根审计并返回 JSON 友好报告。

    ``seed_deltas`` 可用于重放外部记录的 MGP 更新点/轨迹末端；``seed_states``
    还可提供完整的物理网络角/电压状态，使联合求根从已知网络分支开始。两类
    显式种子不足的启动数由固定随机种子补齐。报告中的
    ``used_external_ecritical`` 永远为 ``False``，因为本函数只重算势能，不读取
    外部临界能量。
    """

    if int(max_starts) <= 0:
        raise ValueError("max_starts must be positive")
    if float(root_tol) <= 0 or float(residual_tol) <= 0:
        raise ValueError("root_tol and residual_tol must be positive")

    preset, post = static.preset, static.postfault
    masses = np.asarray(preset.m, dtype=float).reshape(-1)
    sep_delta = _project_coi(np.asarray(post.sep_delta, dtype=float), masses)
    yfull = np.asarray(post.metadata.get("yfull_mod", post.yfull), dtype=complex)
    ngen = int(masses.size)
    nnet = int(yfull.shape[0] - ngen)
    sep_state, sep_ok, sep_residual = spm_energy.solve_spm_network(
        sep_delta, yfull, np.asarray(preset.epu, dtype=float), tol=1e-12,
    )
    if not sep_ok:
        raise RuntimeError(f"SEP network solve failed: residual={sep_residual:g}")

    seeds: list[tuple[np.ndarray, np.ndarray | None, str]] = []
    if seed_deltas is not None:
        for value in seed_deltas:
            arr = np.asarray(value, dtype=float).reshape(-1)
            if arr.size == ngen - 1:
                arr = np.r_[-np.dot(masses[1:], arr) / masses[0], arr]
            if arr.size != ngen or not np.all(np.isfinite(arr)):
                raise ValueError("seed_deltas contains an incompatible or non-finite angle vector")
            seeds.append((_project_coi(arr, masses), None, "explicit_delta"))
    if seed_states is not None:
        for value in seed_states:
            delta, state = _coerce_seed_state(value, ngen, nnet, masses)
            seeds.append((delta, state, "explicit_state"))
    if len(seeds) < int(max_starts):
        generated = _default_seeds(sep_delta, masses, int(max_starts) - len(seeds), random_seed)
        seeds.extend((item, None, "random_delta") for item in generated)
    seeds = seeds[:int(max_starts)]

    records: list[dict] = []
    unique: list[dict] = []
    for index, (seed, network_seed, seed_source) in enumerate(seeds):
        initial = np.r_[seed, sep_state if network_seed is None else network_seed]
        try:
            solution = least_squares(
                lambda u: _joint_residual(u, preset, yfull), initial,
                xtol=float(root_tol), ftol=float(root_tol), gtol=float(root_tol),
                max_nfev=1000,
            )
            raw_delta = _project_coi(solution.x[:ngen], masses)
            delta, wraps = _unwrap_to_reference(raw_delta, sep_delta, masses)
            net = np.asarray(solution.x[ngen:], dtype=float)
            residual_vec = _joint_residual(np.r_[delta, net], preset, yfull)
            network_residual = float(np.linalg.norm(
                spm_energy.spm_network_residual(net, delta, yfull, preset.epu)
            ))
            pe = spm_energy.spm_generator_power(
                delta, net[:nnet], net[nnet:], yfull, np.asarray(preset.epu, dtype=float)
            )
            mismatch = np.asarray(preset.pmpu, dtype=float) - pe
            omega = float(np.sum(mismatch) / np.sum(np.asarray(preset.d, dtype=float)))
            equilibrium_residual = float(np.linalg.norm(mismatch - np.asarray(preset.d, dtype=float) * omega))
            voltage = net[nnet:]
            branch_state, branch_step, branch_ok = _network_branch_to(
                static, delta, sep_delta, sep_state,
            )
            energy = spm_energy.spm_potential_energy(
                preset, post, yfull, sep_delta, sep_state[:nnet], sep_state[nnet:],
                delta, net[:nnet], voltage,
            )
            physical = bool(np.all(np.isfinite(voltage)) and np.all(voltage > 1e-4))
            converged = bool(
                solution.success and np.all(np.isfinite(residual_vec))
                and network_residual < residual_tol
                and equilibrium_residual < residual_tol
                and physical
            )
            record = {
                "start_index": int(index),
                "seed_source": seed_source,
                "solver_success": bool(solution.success),
                "converged": converged,
                "delta_gen": delta.tolist(),
                "theta_net": net[:nnet].tolist(),
                "voltage_net": voltage.tolist(),
                "angle_wraps": wraps,
                "omega_coi": omega,
                "network_residual": network_residual,
                "equilibrium_residual": equilibrium_residual,
                "branch_continuity_error": float(branch_step),
                "branch_continuous": bool(branch_ok),
                "voltage_positive": physical,
                "equilibrium_type": _classify(delta, static),
                "e_critical": float(np.sum(energy)),
                "energy_components": np.asarray(energy, dtype=float).tolist(),
                "exit_reason": "joint root accepted" if converged else "joint root rejected",
            }
            records.append(record)
            if converged and not any(
                np.max(np.abs(np.asarray(record["delta_gen"]) - np.asarray(item["delta_gen"]))) < 1e-6
                for item in unique
            ):
                unique.append(record)
        except Exception as exc:  # noqa: BLE001
            records.append({
                "start_index": int(index), "seed_source": seed_source,
                "solver_success": False, "converged": False,
                "angle_wraps": [], "network_residual": float("nan"),
                "equilibrium_residual": float("nan"), "branch_continuity_error": float("nan"),
                "voltage_positive": False, "branch_continuous": False,
                "equilibrium_type": "UNVERIFIED", "e_critical": float("nan"),
                "energy_components": [], "exit_reason": f"audit failed: {type(exc).__name__}: {exc}",
            })

    return {
        "schema_version": "1.0",
        "audit": "spm_joint_equilibrium_roots",
        "random_seed": int(random_seed),
        "starts_requested": int(max_starts),
        "starts_evaluated": len(records),
        "converged_start_count": int(sum(bool(item.get("converged")) for item in records)),
        "unique_converged_root_count": len(unique),
        "sep_network_residual": float(sep_residual),
        "used_external_ecritical": False,
        "angle_wraps": [item.get("angle_wraps", []) for item in unique],
        "roots": unique,
        "records": records,
        "limitations": [
            "有限初值搜索不是周期状态空间的全根证明",
            "能量由本仓库势能函数重算，不读取外部 E_critical",
            "网络分支连续性只针对 SEP 到候选角的登记路径",
        ],
    }


__all__ = ["audit_spm_joint_roots"]
