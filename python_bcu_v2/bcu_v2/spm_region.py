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


def _branch_id(delta: np.ndarray, theta: np.ndarray, voltage: np.ndarray) -> str:
    payload = np.round(np.r_[delta, theta, voltage], 8).tobytes()
    return hashlib.sha256(payload).hexdigest()[:16]


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

    found: list[SpmEquilibrium] = []
    for delta in candidates:
        if delta.size == preset.ngen - 1:
            delta = np.r_[-np.dot(preset.m[1:], delta) / preset.m[0], delta]
        delta = delta - np.dot(preset.m, delta) / np.sum(preset.m)
        x, ok, network_residual = spm_energy.solve_spm_network(delta, yfull, preset.epu)
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
                found.append(SpmEquilibrium(delta, theta, voltage, eq_type, residual, bid))
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
            stable = bool(result.converged and np.nanmax(np.abs(result.delta_gen)) < 2.0 * np.pi)
            out.append({"clear_time": float(clear), "stable": stable,
                        "status": "APPROXIMATE", "max_algebraic_residual":
                        float(np.nanmax(result.algebraic_residual))})
        except Exception as exc:  # noqa: BLE001
            out.append({"clear_time": float(clear), "stable": None,
                        "status": "BLOCKED", "failure_reason": str(exc)})
    return out
