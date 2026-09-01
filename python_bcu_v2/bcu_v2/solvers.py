# -*- coding: utf-8 -*-
"""P1.1: SciPy 求解器统一层(向后兼容, 可开关).

使用方法:
    - 非线性方程: nlsolve(func, x0, method='scipy'|'newton') 统一接口, scipy 用 optimize.root,
      newton 用 v1 自写牛顿; 可传解析雅可比 jac.
    - SEP: solve_sep_scipy 用 scipy 求平衡点(无副作用, 不改 state).
    - 时域积分: integrate 封装 scipy.solve_ivp(自适应步长/事件, 见 cct.py).
    - benchmark_sep: 同一 SEP 上对比 v1 牛顿 vs scipy(精度/迭代/稳健).

设计: 默认行为向后兼容 v1; scipy 作为可选、更稳更快的后端.
依赖: scipy.
"""

from __future__ import annotations

import time
from typing import Callable, Optional, Tuple

import numpy as np


def nlsolve(func: Callable[[np.ndarray], np.ndarray], x0: np.ndarray,
            method: str = "scipy", jac=None, tol: float = 1e-10,
            max_iter: int = 200) -> Tuple[np.ndarray, bool, float]:
    """统一非线性求解: 返回 (解, 是否收敛, 残差范数).

    使用方法: method='scipy' 用 optimize.root(hybr, 可选 jac); method='newton' 用 v1 牛顿.
    """

    if method == "scipy":
        from scipy.optimize import root
        sol = root(func, np.asarray(x0, dtype=float), jac=jac, method="hybr", tol=tol)
        r = float(np.linalg.norm(np.asarray(func(sol.x), dtype=float)))
        return np.asarray(sol.x, dtype=float), bool(sol.success and r < max(1e3 * tol, 1e-6)), r
    elif method == "newton":
        from bcu_3m9b.numerics import newton_solve
        x, ok, _, r = newton_solve(func, x0, tol=tol, max_iter=max_iter)
        return x, bool(ok), float(r)
    raise ValueError(f"unknown method {method!r}")


def solve_sep_scipy(preset, state, basevalue, delta0: Optional[np.ndarray] = None,
                    omega0: float = 0.0, tol: float = 1e-10) -> Tuple[np.ndarray, float, np.ndarray, bool]:
    """用 scipy 求 SEP(无副作用, 不改写 state).

    使用方法: 返回 (COI 角度, pu 速度, 功率残差, 是否收敛), 与 v1.solve_sep 数值等价.
    """

    from bcu_3m9b.equilibrium import sep_residual, normalize_coi

    n = preset.ngen
    if delta0 is None:
        delta0 = np.zeros(n)
    delta0 = normalize_coi(delta0, preset.m)
    z0 = np.r_[delta0[: n - 1] - delta0[-1], omega0]
    x, ok, _ = nlsolve(lambda z: sep_residual(z, preset, state, basevalue), z0, tol=tol)
    delta = normalize_coi(np.r_[x[: n - 1], 0.0], preset.m)
    omegapu = 1.0 + x[-1] / basevalue.omega_b
    perr = sep_residual(x, preset, state, basevalue)
    return delta, omegapu, perr, bool(ok and np.linalg.norm(perr) <= max(10 * tol, 1e-7))


def integrate(rhs: Callable, x0: np.ndarray, t_span, t_eval=None,
              method: str = "RK45", rtol: float = 1e-8, atol: float = 1e-10,
              events=None):
    """封装 scipy.solve_ivp(自适应步长, 可事件检测).

    使用方法: 传入 rhs(t,x), 初值, 时间区间; 返回 solve_ivp 的解对象.
    """

    from scipy.integrate import solve_ivp
    return solve_ivp(rhs, t_span, np.asarray(x0, dtype=float), t_eval=t_eval,
                     method=method, rtol=rtol, atol=atol, events=events)


def benchmark_sep(static, bad_offset: float = 1.0) -> dict:
    """对比 v1 牛顿 vs scipy 求同一 SEP: 精度一致、以及从坏初值的稳健性.

    使用方法: 返回 dict, 含两法解的差、各自是否从坏初值收敛、耗时.
    """

    from bcu_3m9b.equilibrium import solve_sep

    preset, post, base = static.preset, static.postfault, static.basevalue
    good = np.zeros(preset.ngen)

    t0 = time.perf_counter()
    d_v1, _, _, ok_v1, _ = solve_sep(preset, post, base, good.copy(), 0.0, tol=1e-10)
    t_v1 = time.perf_counter() - t0
    # v1 会改写 post.SEP, 复原以免污染后续.
    post.sep_delta = d_v1

    t0 = time.perf_counter()
    d_sp, _, _, ok_sp = solve_sep_scipy(preset, post, base, good.copy(), 0.0, tol=1e-10)
    t_sp = time.perf_counter() - t0

    match = float(np.max(np.abs(d_v1 - d_sp)))

    # 稳健性: 从坏初值出发.
    bad = good + bad_offset
    _, _, _, ok_v1_bad, _ = solve_sep(preset, post, base, bad.copy(), 0.0, tol=1e-10)
    post.sep_delta = d_v1  # 再复原
    _, _, _, ok_sp_bad = solve_sep_scipy(preset, post, base, bad.copy(), 0.0, tol=1e-10)

    return {"match": match, "ok_v1": bool(ok_v1), "ok_scipy": bool(ok_sp),
            "ok_v1_badguess": bool(ok_v1_bad), "ok_scipy_badguess": bool(ok_sp_bad),
            "time_v1": t_v1, "time_scipy": t_sp}
