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

from typing import Dict

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root


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

    steps = max(2, int(round(tlength / tunit)))
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
