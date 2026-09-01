"""SPM 平衡点与稳定域抽样接口。

平衡点角度候选沿用 v2 reduced 梯度搜索，再对每个候选求完整 SPM 网络状态；输出包含
残差、类型和稳定分支 ID，避免只比较点数。区域抽样在 MATLAB 参考可用前保持
``APPROXIMATE``/``UNVERIFIED`` 标记。
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib

import numpy as np

from . import spm_energy


@dataclass
class SpmEquilibrium:
    delta_gen: np.ndarray
    theta_net: np.ndarray
    voltage_net: np.ndarray
    equilibrium_type: str
    residual_norm: float
    branch_id: str
    continuity_error: float = 0.0


def _branch_id(delta: np.ndarray, theta: np.ndarray, voltage: np.ndarray) -> str:
    payload = np.round(np.r_[delta, theta, voltage], 8).tobytes()
    return hashlib.sha256(payload).hexdigest()[:16]


def _solve_network_continuous(target_delta: np.ndarray, sep_delta: np.ndarray,
                              yfull: np.ndarray, epu: np.ndarray,
                              sep_state: np.ndarray, *, segments: int = 50):
    """沿 SEP 到目标角的直线连续求网络代数分支。

    SPM 网络方程在同一发电机角下可能有多个数学根。每个候选平衡点都
    从已经验证的 SEP 网络状态 warm-start，并逐段校正；冷启动零/一初值
    不能作为区域比较的分支选择器。返回 ``(state, max_step, residual)``，
    失败时返回 ``(None, inf, inf)``。
    """

    if segments < 1:
        raise ValueError("segments must be positive")
    current = np.asarray(sep_state, dtype=float).copy()
    max_step = 0.0
    previous = current.copy()
    for alpha in np.linspace(1.0 / segments, 1.0, segments):
        delta = (1.0 - alpha) * np.asarray(sep_delta, dtype=float) + alpha * np.asarray(target_delta, dtype=float)
        state, ok, residual = spm_energy.solve_spm_network(
            delta, yfull, epu, guess=current, tol=1e-11
        )
        if not ok or not np.all(np.isfinite(state)):
            return None, float("inf"), float(residual)
        max_step = max(max_step, float(np.linalg.norm(state - previous)))
        previous = np.asarray(state, dtype=float)
        current = previous
    return current, max_step, float(np.linalg.norm(
        spm_energy.spm_network_residual(
            current, np.asarray(target_delta, dtype=float), yfull, epu
        )
    ))


def enumerate_spm_equilibria(static, *, grid_points: int = 21,
                             include_sep: bool = True) -> list[SpmEquilibrium]:
    """枚举 SPM 网络平衡点并按角度/网络状态去重。"""

    from bcu_3m9b.experiments import find_reduced_equilibria

    preset, post = static.preset, static.postfault
    yfull = np.asarray(post.metadata.get("yfull_mod", post.yfull), dtype=complex)
    candidates = []
    if include_sep:
        candidates.append(np.asarray(post.sep_delta, dtype=float))
    for item in find_reduced_equilibria(post, preset, grid_points=grid_points):
        candidates.append(np.asarray(item["xep"], dtype=float))

    sep_state, sep_ok, sep_network_residual = spm_energy.solve_spm_network(
        np.asarray(post.sep_delta, dtype=float), yfull, preset.epu
    )
    if not sep_ok:
        return []

    found: list[SpmEquilibrium] = []
    for delta in candidates:
        if delta.size == preset.ngen - 1:
            delta = np.r_[-np.dot(preset.m[1:], delta) / preset.m[0], delta]
        delta = delta - np.dot(preset.m, delta) / np.sum(preset.m)
        if np.max(np.abs(delta - np.asarray(post.sep_delta, dtype=float))) < 1e-8:
            x, ok, network_residual = sep_state, True, sep_network_residual
            continuity_error = 0.0
        else:
            x, continuity_error, network_residual = _solve_network_continuous(
                delta, np.asarray(post.sep_delta, dtype=float), yfull, preset.epu,
                sep_state,
            )
            ok = x is not None and np.isfinite(network_residual) and network_residual < 1e-6
        if not ok:
            continue
        theta, voltage = x[:preset.nbus - preset.ngen], x[preset.nbus - preset.ngen:]
        pe = spm_energy.spm_generator_power(delta, theta, voltage, yfull, preset.epu)
        mismatch = preset.pmpu - pe
        pcoi = np.sum(mismatch)
        projected = mismatch - preset.m / np.sum(preset.m) * pcoi
        residual = float(np.linalg.norm(np.r_[projected, network_residual]))
        if residual < 1e-6:
            try:
                from .cuep import _reduced_jacobian_eig
                eig = _reduced_jacobian_eig(delta, post.yred, preset)
                npos = int(np.sum(eig > 1e-6))
                eq_type = "SEP" if npos == 0 else ("type-1" if npos == 1 else f"type-{npos}")
            except Exception:
                eq_type = "unknown"
            bid = _branch_id(delta, theta, voltage)
            if not any(np.max(np.abs(delta - e.delta_gen)) < 1e-6 for e in found):
                found.append(SpmEquilibrium(delta, theta, voltage, eq_type, residual, bid,
                                            float(continuity_error)))
    found.sort(key=lambda item: (item.equilibrium_type != "SEP", item.branch_id))
    return found


def sample_spm_region(static, *, grid_points: int = 9,
                      clear_times: np.ndarray | None = None) -> list[dict]:
    """以有限清除时间抽样稳定性；结果明确标为近似。"""

    from .spm_dae import simulate_spm_trajectory

    times = np.asarray(clear_times if clear_times is not None else
                       np.linspace(0.02, 0.2, grid_points), dtype=float)
    out = []
    for clear in times:
        try:
            result = simulate_spm_trajectory(static, clear_time=float(clear),
                                             postfault_time=0.3, tunit=0.005)
            finite_residual = np.all(np.isfinite(result.algebraic_residual))
            stable = bool(result.converged and finite_residual and
                          np.max(np.abs(result.delta_gen)) < 2.0 * np.pi)
            classification = "stable" if stable else "unstable"
            out.append({"clear_time": float(clear), "stable": stable,
                        "classification": classification,
                        "status": "APPROXIMATE", "max_algebraic_residual":
                        float(np.max(result.algebraic_residual)) if finite_residual else float("nan")})
        except Exception as exc:  # noqa: BLE001
            out.append({"clear_time": float(clear), "stable": None,
                        "status": "BLOCKED", "failure_reason": str(exc)})
    return out
