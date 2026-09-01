"""基础数值工具使用说明。

使用方法：
    本模块只依赖 NumPy,不需要 SciPy。求解非线性方程用 ``newton_solve``;
    需要雅可比时用 ``numerical_jacobian``;做一步定步长积分用 ``rk4_step``;
    角度归一化用 ``wrap_angle``。

设计说明：
    为让平台在无 SciPy 的纯 Python 环境也能启动,这里自带有限差分雅可比、
    阻尼牛顿法和四阶 RK4。
"""

from __future__ import annotations

from typing import Callable, Tuple

import numpy as np


def numerical_jacobian(func: Callable[[np.ndarray], np.ndarray], x: np.ndarray,
                       step: float = 1e-6) -> np.ndarray:
    """用中心差分近似向量函数的雅可比矩阵。

    使用方法：
        传入向量函数 ``func`` 和当前点 ``x``,返回 f 对 x 的雅可比矩阵。
    参数：
        func：接受长度 n 向量、返回长度 m 向量的函数。
        x：求导所在的点。
        step：相对差分步长。
    返回：
        形状 (m, n) 的雅可比矩阵。
    步骤：
        对每个分量正负各扰动一次,用 (f(x+h)-f(x-h))/(2h) 填入对应列。
    """

    x = np.asarray(x, dtype=float)
    f0 = np.asarray(func(x), dtype=float)
    jac = np.zeros((f0.size, x.size), dtype=float)
    for i in range(x.size):
        # 步长按分量幅值缩放,避免大分量处相对精度不足。
        h = step * max(1.0, abs(x[i]))
        xp = x.copy()
        xm = x.copy()
        xp[i] += h
        xm[i] -= h
        jac[:, i] = (np.asarray(func(xp)) - np.asarray(func(xm))) / (2.0 * h)
    return jac


def newton_solve(func: Callable[[np.ndarray], np.ndarray], x0: np.ndarray,
                 tol: float = 1e-10, max_iter: int = 100,
                 jac_step: float = 1e-6, damping: bool = True,
                 name: str = "nonlinear system") -> Tuple[np.ndarray, bool, int, float]:
    """阻尼牛顿法求解 f(x)=0。

    使用方法：
        传入残差函数 ``func`` 和初值 ``x0``,返回 (解, 是否收敛, 迭代次数, 残差范数)。
        需要更严/更松时调 ``tol``;病态雅可比会自动退回最小二乘步。
    参数：
        func：残差向量函数。
        x0：迭代初值。
        tol：收敛阈值（残差二范数）。
        max_iter：最大迭代步数。
        jac_step：数值雅可比的差分步长。
        damping：是否启用回溯阻尼（保证残差单调下降）。
        name：诊断标签,便于区分不同方程。
    返回：
        (x, success, iterations, residual_norm)。
    步骤：
        每步先算雅可比与牛顿步,启用阻尼时按 0.5 递减步长直到残差下降,
        否则直接整步;步长与残差同时达标即判定收敛。
    """

    x = np.asarray(x0, dtype=float).reshape(-1).copy()
    f = np.asarray(func(x), dtype=float).reshape(-1)
    norm_f = float(np.linalg.norm(f, ord=2))
    for iteration in range(1, max_iter + 1):
        if not np.isfinite(norm_f):
            return x, False, iteration, norm_f
        if norm_f <= tol:
            return x, True, iteration - 1, norm_f
        jac = numerical_jacobian(func, x, jac_step)
        try:
            step = np.linalg.solve(jac, -f)
        except np.linalg.LinAlgError:
            # 雅可比奇异时退回最小二乘解,保证迭代能继续推进。
            step = np.linalg.lstsq(jac, -f, rcond=None)[0]
        alpha = 1.0
        accepted = False
        if damping:
            # 回溯线搜索：不断折半步长,接受第一个让残差下降的试探点。
            for _ in range(12):
                trial = x + alpha * step
                trial_f = np.asarray(func(trial), dtype=float).reshape(-1)
                trial_norm = float(np.linalg.norm(trial_f, ord=2))
                if np.isfinite(trial_norm) and trial_norm < norm_f:
                    accepted = True
                    x, f, norm_f = trial, trial_f, trial_norm
                    break
                alpha *= 0.5
        if not accepted:
            x = x + step
            f = np.asarray(func(x), dtype=float).reshape(-1)
            norm_f = float(np.linalg.norm(f, ord=2))
        if np.linalg.norm(step) * alpha <= tol and norm_f <= 10 * tol:
            return x, True, iteration, norm_f
    return x, norm_f <= tol, max_iter, norm_f


def rk4_step(func: Callable[[np.ndarray], np.ndarray], x: np.ndarray,
             step: float) -> np.ndarray:
    """经典四阶 Runge--Kutta 单步推进。

    使用方法：
        传入右端函数 ``func``、当前状态 ``x`` 和步长 ``step``,返回下一步状态。
    步骤：
        依次求 k1..k4 四个斜率,按 (k1+2k2+2k3+k4)/6 加权得到增量。
    """

    k1 = np.asarray(func(x), dtype=float)
    k2 = np.asarray(func(x + 0.5 * step * k1), dtype=float)
    k3 = np.asarray(func(x + 0.5 * step * k2), dtype=float)
    k4 = np.asarray(func(x + step * k3), dtype=float)
    return x + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0


def wrap_angle(angle: np.ndarray) -> np.ndarray:
    """把角度归一化到 [-pi, pi)。

    使用方法：
        传入任意角度（rad）,返回等价的主值角度,便于比较和作图。
    """

    return (np.asarray(angle) + np.pi) % (2.0 * np.pi) - np.pi
