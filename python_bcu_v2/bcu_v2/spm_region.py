"""SPM 平衡点与稳定域抽样接口。

平衡点角度候选沿用 v2 reduced 梯度搜索，再对每个候选求完整 SPM 网络状态；输出包含
残差、类型和稳定分支 ID，避免只比较点数。固定稳定流形检查点使用约束切空间积分，
可与 MATLAB 原生参考逐变量对照；更广泛的清除时间区域抽样仍明确标为 ``APPROXIMATE``。
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib

import numpy as np
from scipy.integrate import solve_ivp

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


def branch_continuity(states: np.ndarray, *, tolerance: float = 1.0) -> tuple[bool, float]:
    """检查一组连续网络状态是否发生不可接受的分支跳跃。

    返回 ``(ok, max_step)``；空/非有限状态直接失败。该检查同时用于
    warm-start 诊断和变异测试，不能通过关闭检查来掩盖分支切换。
    """

    values = np.asarray(states, dtype=float)
    if values.ndim != 2 or values.shape[0] < 2 or not np.all(np.isfinite(values)):
        return False, float("inf")
    max_step = float(np.max(np.linalg.norm(np.diff(values, axis=0), axis=1)))
    return bool(max_step <= float(tolerance)), max_step


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


def trace_spm_stable_manifold(static, equilibrium: SpmEquilibrium,
                              perturb_vector: np.ndarray, sign: float, *,
                              perturb: float = 1e-2,
                              sample_times: np.ndarray | None = None,
                              rtol: float = 1e-10,
                              atol: float = 1e-12) -> dict:
    """Trace one stable-manifold branch on the algebraic network manifold.

    MATLAB's ``Statable_Region_SPM`` integrates a stiff DAE backwards from a
    type-1 UEP.  This implementation keeps the same two reduced generator
    coordinates but solves the SPM network algebraic equations continuously at
    every RHS evaluation.  No zero/one cold-start state is used after the
    initial physical equilibrium branch has been established.

    The returned dictionary is deliberately explicit about convergence and
    residuals so callers can keep a partial/blocked comparison separate from
    the equilibrium-point comparison.
    """

    if equilibrium.equilibrium_type != "type-1":
        raise ValueError("stable-manifold tracing requires a type-1 equilibrium")
    if perturb <= 0 or not np.isfinite(perturb):
        raise ValueError("perturb must be positive and finite")
    times = np.asarray(sample_times if sample_times is not None
                       else np.array([0.0, 0.25, 0.5, 0.75, 1.0]), dtype=float)
    if times.ndim != 1 or times.size < 2 or not np.all(np.isfinite(times)):
        raise ValueError("sample_times must be a finite one-dimensional array")
    if np.any(np.diff(times) <= 0):
        raise ValueError("sample_times must be strictly increasing")
    vector = np.asarray(perturb_vector, dtype=float).reshape(-1)
    if vector.size != 2 or not np.all(np.isfinite(vector)):
        raise ValueError("perturb_vector must contain two finite values")
    norm = float(np.linalg.norm(vector))
    if norm == 0:
        raise ValueError("perturb_vector must be nonzero")
    vector = vector / norm
    sign_value = float(np.sign(sign))
    if sign_value == 0:
        raise ValueError("sign must be nonzero")

    preset = static.preset
    post = static.postfault
    yfull = np.asarray(post.metadata.get("yfull_mod", post.yfull), dtype=complex)
    epu = np.asarray(preset.epu, dtype=float)
    ngen = int(preset.ngen)
    nnet = int(yfull.shape[0] - ngen)
    m = np.asarray(preset.m, dtype=float)
    sep_state, sep_ok, sep_residual = spm_energy.solve_spm_network(
        np.asarray(equilibrium.delta_gen, dtype=float), yfull, epu,
        guess=np.r_[equilibrium.theta_net, equilibrium.voltage_net],
        tol=1e-12,
    )
    if not sep_ok:
        return {"time": times, "converged": False,
                "failure_reason": f"UEP network warm-start failed: residual={sep_residual:g}"}

    x0 = np.asarray(equilibrium.delta_gen, dtype=float)[1:] + sign_value * perturb * vector
    cache = {"z": sep_state.copy()}

    def make_delta(x: np.ndarray) -> np.ndarray:
        d23 = np.asarray(x, dtype=float).reshape(-1)
        d1 = -float(np.dot(m[1:], d23) / m[0])
        d = np.r_[d1, d23]
        return d - np.dot(m, d) / np.sum(m)

    def network_residual(z: np.ndarray, delta: np.ndarray) -> np.ndarray:
        return spm_energy.spm_network_residual(z, delta, yfull, epu)

    def finite_jacobian(fun, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        f0 = np.asarray(fun(x), dtype=float)
        jac = np.empty((f0.size, x.size), dtype=float)
        for j in range(x.size):
            step = 1e-6 * max(1.0, abs(float(x[j])))
            xp = x.copy(); xm = x.copy()
            xp[j] += step; xm[j] -= step
            jac[:, j] = (np.asarray(fun(xp)) - np.asarray(fun(xm))) / (2.0 * step)
        return jac

    # Integrate the differential-algebraic tangent explicitly after
    # differentiating g(delta,z)=0.  This is the index-1 constrained ODE
    # represented by MATLAB's tiny algebraic mass entries and avoids Newton
    # branch jumps during implicit solver stage evaluations.
    initial_delta = make_delta(x0)
    z_initial, initial_ok, initial_residual = spm_energy.solve_spm_network(
        initial_delta, yfull, epu, guess=sep_state, tol=1e-11
    )
    if not initial_ok:
        return {"time": times, "converged": False,
                "failure_reason": f"perturbed UEP network solve failed: residual={initial_residual:g}"}
    z_initial = np.asarray(z_initial, dtype=float)

    def rhs(_t: float, state_vec: np.ndarray) -> np.ndarray:
        x = np.asarray(state_vec[:2], dtype=float)
        z = np.asarray(state_vec[2:], dtype=float)
        delta = make_delta(x)
        pe = spm_energy.spm_generator_power(delta, z[:nnet], z[nnet:], yfull, epu)
        mismatch = np.asarray(preset.pmpu, dtype=float) - pe
        projected = mismatch - m / np.sum(m) * np.sum(mismatch)
        # f_reducedstate_SPM_backward uses the negative projected generator
        # acceleration for the two independent COI coordinates.
        delta_dot = -projected[1:] / m[1:]
        gz = finite_jacobian(lambda zz: network_residual(zz, delta), z)
        gd = finite_jacobian(lambda dd: network_residual(z, make_delta(dd)), x)
        try:
            z_dot = np.linalg.solve(gz, -(gd @ delta_dot))
        except np.linalg.LinAlgError as exc:
            raise RuntimeError("stable-manifold constraint Jacobian is singular") from exc
        return np.r_[delta_dot, z_dot]

    try:
        sol = solve_ivp(rhs, (float(times[0]), float(times[-1])),
                        np.r_[x0, z_initial],
                        method="RK45", t_eval=times, rtol=rtol, atol=atol,
                        # Keep continuation increments small enough that the
                        # constrained tangent remains on the positive-voltage
                        # branch between output points.
                        max_step=min(5e-3, float(np.min(np.diff(times)))) )
    except Exception as exc:  # noqa: BLE001
        return {"time": times, "converged": False, "failure_reason": str(exc)}
    if not sol.success or sol.y.shape[1] != times.size:
        return {"time": times, "converged": False,
                "failure_reason": str(sol.message)}

    delta_out = np.vstack([make_delta(sol.y[:2, k]) for k in range(times.size)])
    theta_out = np.empty((times.size, nnet), dtype=float)
    voltage_out = np.empty((times.size, nnet), dtype=float)
    residuals = np.empty(times.size, dtype=float)
    states = []
    for k, delta in enumerate(delta_out):
        z = np.asarray(sol.y[2:, k], dtype=float)
        states.append(z.copy())
        theta_out[k] = z[:nnet]
        voltage_out[k] = z[nnet:]
        residuals[k] = float(np.linalg.norm(network_residual(z, delta)))
    continuity = float(max(np.linalg.norm(np.diff(np.asarray(states), axis=0), axis=1), default=0.0))
    return {"time": times, "delta_gen": delta_out,
            "theta_net": theta_out, "voltage_net": voltage_out,
            "residual_norm": residuals, "branch_continuity_error": continuity,
            "converged": bool(np.all(np.isfinite(residuals))), "sign": sign_value}


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
