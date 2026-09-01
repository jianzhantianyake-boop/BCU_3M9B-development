# -*- coding: utf-8 -*-
"""P1.3: 严格 DAE 级的结构保持模型(SPM), 根治代数解脆弱.

思路(约束流形上的降阶 DAE):
    结构保持模型是 index-1 DAE: 发电机 (δg, ωg) 为微分变量, 负荷母线 (δL, V) 为代数变量,
    受功率平衡约束 g(δg, z)=0 约束. 这里把它当"约束流形上的降阶 ODE"求解:
      - 用 scipy.solve_ivp 自适应积分发电机状态(method='RK45' 或刚性 'Radau');
      - 每次 RHS 求值时用**连续法**(上一步解热启动的 scipy.root 校正)把 g=0 解到机器精度,
        使代数约束在每个求值点都严格满足. 相比 v1"每步一次冷启动牛顿", 连续法的初值始终贴近
        真解, 从而根治偶发不收敛.

对比 v1:
    v1.simulate_spm: 固定步 + 每步冷启动/上一步热启动牛顿, 偶发首步不收敛即抛异常.
    本模块: 自适应步 + 连续法热启动 + scipy.root + 失败回退(default guess / 同伦), 更稳更准.

依赖: scipy.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root


@dataclass
class SpmTrajectoryResult:
    """统一的 SPM 四阶段轨迹结果。"""

    time: np.ndarray
    delta_gen: np.ndarray
    omega_gen: np.ndarray
    theta_net: np.ndarray
    voltage_net: np.ndarray
    algebraic_residual: np.ndarray
    phase_labels: np.ndarray
    converged: bool


def _algebraic_context(state, preset):
    """内部: 取 SPM 代数方程所需的重排导纳 / 负荷 / 维度."""

    from bcu_3m9b.spm import _load_power

    yorg = np.asarray(state.metadata.get("yorg_mod", state.yfull), dtype=complex)
    transform = np.asarray(state.metadata.get("transform"), dtype=int)
    n = preset.ngen
    load_pq = _load_power(preset, transform[n:])
    nload = yorg.shape[0] - n
    return yorg, load_pq, n, nload


def _make_solver(yorg, load_pq, ngen, nload, tol=1e-11):
    """内部: 返回带连续法热启动 + 回退的代数求解闭包 solve(δg)->z."""

    from bcu_3m9b.spm import algebraic_residual

    default_guess = np.r_[np.zeros(nload), np.ones(nload)]
    cache = {"z": default_guess.copy()}

    def solve(delta_gen: np.ndarray) -> np.ndarray:
        dg = np.asarray(delta_gen, dtype=float)

        def resid(z):
            return algebraic_residual(z, dg, yorg, load_pq, ngen)

        # 校正: 先用上一步解(连续法热启动).
        sol = root(resid, cache["z"], method="hybr", tol=tol)
        z, r = sol.x, float(np.linalg.norm(resid(sol.x)))
        if r > 1e-7:  # 回退1: 默认初值
            s2 = root(resid, default_guess, method="hybr", tol=tol)
            if np.linalg.norm(resid(s2.x)) < r:
                z, r = s2.x, float(np.linalg.norm(resid(s2.x)))
        if r > 1e-7:  # 回退2: lm(阻尼最小二乘)
            s3 = root(resid, cache["z"], method="lm")
            if np.linalg.norm(resid(s3.x)) < r:
                z, r = s3.x, float(np.linalg.norm(resid(s3.x)))
        cache["z"] = z
        return z

    def reset(z=None):
        cache["z"] = default_guess.copy() if z is None else np.asarray(z, dtype=float)

    solve.reset = reset  # type: ignore[attr-defined]
    return solve


def simulate_spm_dae(tlength: float, tunit: float, state, preset, basevalue,
                     delta0: np.ndarray, omega0: np.ndarray,
                     method: str = "RK45", rtol: float = 1e-8, atol: float = 1e-10) -> Dict:
    """严格 DAE 级 SPM 仿真(约束流形降阶 ODE + 连续法).

    使用方法:
        传入时长/输出步长/网络工况/参数/基值/初值; method 可选 'RK45'(默认, 快)或
        'Radau'/'BDF'(刚性). 返回 {time, delta, omega, delta_coi, algebraic, success, method}.
    """

    from bcu_3m9b.spm import spm_generator_rhs

    yorg, load_pq, n, nload = _algebraic_context(state, preset)
    solve_alg = _make_solver(yorg, load_pq, n, nload)

    # 一致初始化: 在初始发电机角上解代数, 作为流形起点.
    solve_alg.reset()
    z0 = solve_alg(np.asarray(delta0, dtype=float))

    def rhs(t, x):
        dg, om = x[:n], x[n:]
        z = solve_alg(dg)  # 连续法: 约束在每个求值点严格满足
        return spm_generator_rhs(dg, om, state, preset, basevalue, z)

    # ``tunit`` is a requested output spacing, so include both endpoints:
    # N intervals require N+1 samples.  The previous N-sample grid drifted
    # every checkpoint by one interval over a long fault trajectory.
    steps = max(2, int(round(tlength / tunit)) + 1)
    t_eval = np.linspace(0.0, tlength, steps)
    sol = solve_ivp(rhs, [0.0, tlength], np.r_[delta0, omega0], method=method,
                    t_eval=t_eval, rtol=rtol, atol=atol, dense_output=False)

    delta = sol.y[:n].T
    omega = sol.y[n:].T
    # 沿输出时刻重建代数量(顺序推进 -> 连续法跟踪).
    solve_alg.reset(z0)
    algebraic = np.array([solve_alg(delta[k]) for k in range(sol.t.size)])
    msum = np.sum(preset.m)
    delta_coi = delta - (delta @ preset.m / msum)[:, None]
    return {"time": sol.t, "delta": delta, "omega": omega, "algebraic": algebraic,
            "delta_coi": delta_coi, "success": bool(sol.success), "method": method}


def simulate_spm_trajectory(static, *, clear_time: float = 0.2,
                            postfault_time: float = 0.3, tunit: float = 1e-3,
                            method: str = "RK45") -> SpmTrajectoryResult:
    """运行 prefault/fault-on/clearing/postfault recovery 四阶段 SPM 轨迹。

    清除时刻沿用故障段最后一个状态；网络代数状态在每个输出点重新校正并报告残差。
    ``RK45`` 与 ``Radau`` 可用同一接口重复运行比较。
    """

    if clear_time <= 0 or postfault_time <= 0 or tunit <= 0:
        raise ValueError("clear_time, postfault_time and tunit must be positive")
    preset, base = static.preset, static.basevalue
    delta0 = np.asarray(static.prefault.sep_delta, dtype=float)
    omega0 = np.full(preset.ngen, float(static.prefault.sep_omegapu) * base.omega_b)
    fault = simulate_spm_dae(clear_time, tunit, static.fault, preset, base,
                              delta0, omega0, method=method)
    delta_clear = fault["delta"][-1]
    omega_clear = fault["omega"][-1]
    post = simulate_spm_dae(postfault_time, tunit, static.postfault, preset, base,
                             delta_clear, omega_clear, method=method)

    # 保留 t=0 的 prefault 标签，故障段从第二个点开始，清除点独立标注。
    f_slice = slice(1, None)
    p_slice = slice(1, None)
    time = np.r_[0.0, np.asarray(fault["time"])[f_slice],
                 clear_time + np.asarray(post["time"])[p_slice]]
    delta = np.vstack([fault["delta"][0], fault["delta"][f_slice], post["delta"][p_slice]])
    omega = np.vstack([fault["omega"][0], fault["omega"][f_slice], post["omega"][p_slice]])
    # fault-on 网络可能移除故障母线，代数量维度小于 postfault；统一输出宽度时
    # 用 NaN 标记“该阶段不存在的节点”，不以 0 冒充可验证状态。
    raw_alg = [np.atleast_2d(fault["algebraic"][0]), fault["algebraic"][f_slice], post["algebraic"][p_slice]]
    nload_parts = [int(a.shape[1] // 2) for a in raw_alg]
    max_nload = max(nload_parts)
    max_alg = 2 * max_nload
    algebraic_parts = []
    theta_parts = []
    voltage_parts = []
    for part, nload_part in zip(raw_alg, nload_parts):
        padded = np.full((part.shape[0], max_alg), np.nan, dtype=float)
        padded[:, :nload_part] = part[:, :nload_part]
        padded[:, max_nload:max_nload + nload_part] = part[:, nload_part:]
        algebraic_parts.append(padded)
        theta_pad = np.full((part.shape[0], max_nload), np.nan, dtype=float)
        voltage_pad = np.full((part.shape[0], max_nload), np.nan, dtype=float)
        theta_pad[:, :nload_part] = part[:, :nload_part]
        voltage_pad[:, :nload_part] = part[:, nload_part:]
        theta_parts.append(theta_pad)
        voltage_parts.append(voltage_pad)
    algebraic = np.vstack(algebraic_parts)
    n_fault = max(0, fault["time"].size - 1)
    n_post = max(0, post["time"].size - 1)
    labels = np.asarray(["prefault"] + ["fault-on"] * max(0, n_fault - 1) +
                        (["clearing"] if n_fault else []) +
                        ["postfault recovery"] * n_post, dtype=object)
    yorg_f, load_f, n, nload_f = _algebraic_context(static.fault, preset)
    yorg_p, load_p, _, nload_p = _algebraic_context(static.postfault, preset)
    residual = []
    for k in range(delta.shape[0]):
        if k <= n_fault:
            z = np.asarray(fault["algebraic"][k], dtype=float)
            value = float(np.linalg.norm(__import__("bcu_3m9b.spm", fromlist=["algebraic_residual"]).algebraic_residual(
                z, delta[k], yorg_f, load_f, n)))
            residual.append(value if value < 1e-7 else np.nan)
        else:
            post_index = k - n_fault
            raw_post = post["algebraic"][post_index]
            z = np.asarray(raw_post, dtype=float)
            value = float(np.linalg.norm(__import__("bcu_3m9b.spm", fromlist=["algebraic_residual"]).algebraic_residual(
                z, delta[k], yorg_p, load_p, n)))
            residual.append(value if value < 1e-7 else np.nan)
    theta_net = np.vstack(theta_parts)
    voltage_net = np.vstack(voltage_parts)
    return SpmTrajectoryResult(time=time, delta_gen=delta, omega_gen=omega,
                               theta_net=theta_net, voltage_net=voltage_net,
                               algebraic_residual=np.asarray(residual),
                               phase_labels=labels,
                               converged=bool(fault["success"] and post["success"] and
                                              np.all(np.isfinite(residual))))


def select_spm_checkpoints(result: SpmTrajectoryResult, *, clear_time: float,
                           postfault_time: float, tunit: float) -> list[dict]:
    """Return fixed, nearest-output SPM checkpoints for cross-validation.

    No interpolation or zero padding is performed.  If a requested time is
    outside the simulated horizon, the record is marked ``available=False``
    instead of silently reusing the final state.  The pre-clearing checkpoint
    is the last output strictly before ``clear_time``; the clearing checkpoint
    must be exactly at the fault trajectory boundary for a valid comparison.
    """
    if clear_time <= 0 or postfault_time <= 0 or tunit <= 0:
        raise ValueError("clear_time, postfault_time and tunit must be positive")
    requested = [
        ("t0", 0.0),
        ("pre-clearing", max(0.0, clear_time - tunit)),
        ("clearing", clear_time),
        ("post-clearing-10ms", clear_time + 0.01),
        ("post-clearing-50ms", clear_time + 0.05),
        ("post-clearing-100ms", clear_time + 0.10),
        ("final", clear_time + postfault_time),
    ]
    times = np.asarray(result.time, dtype=float)
    out: list[dict] = []
    for label, target in requested:
        if times.size == 0 or target < times[0] - 1e-12 or target > times[-1] + 1e-12:
            out.append({"label": label, "requested_time": float(target), "available": False})
            continue
        idx = int(np.argmin(np.abs(times - target)))
        actual = float(times[idx])
        record = {
            "label": label,
            "requested_time": float(target),
            "actual_time": actual,
            "index": idx,
            "available": True,
            "phase": str(result.phase_labels[idx]),
            "delta_gen": np.asarray(result.delta_gen[idx], dtype=float).copy(),
            "omega_gen": np.asarray(result.omega_gen[idx], dtype=float).copy(),
            "theta_net": np.asarray(result.theta_net[idx], dtype=float).copy(),
            "voltage_net": np.asarray(result.voltage_net[idx], dtype=float).copy(),
            "algebraic_residual": float(result.algebraic_residual[idx]),
        }
        out.append(record)
    return out
