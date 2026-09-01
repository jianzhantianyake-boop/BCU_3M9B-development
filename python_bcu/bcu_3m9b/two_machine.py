"""两机模型教学实现使用说明。

使用方法：
    先用 ``TwoMachineParameters`` 显式装入网络和机械参数；算两机电磁功率用
    ``electrical_power``；完整/故障/约简摆动方程分别用 ``f_2m``、``f_2m_fault``、
    ``f_2m_reduce``；扫描并分类平衡点用 ``equilibria``；积分轨迹用
    ``simulate_two_machine``。

对应关系：
    对应 MATLAB ``f_2m.m``、``f_2m_fault.m`` 和 ``f_2m_reduce.m``，参数不再从 base
    workspace 读取，而是通过参数对象显式传入。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Tuple

import numpy as np

from .numerics import newton_solve, numerical_jacobian, rk4_step


@dataclass
class TwoMachineParameters:
    """两机等值网络与机械参数容器。

    使用方法：
        构造时按顺序传入电导/电纳、内电势、机械功率、惯量和阻尼；故障电导
        ``fault_g1``、``fault_g2`` 缺省时故障期沿用正常电导。
    """

    g1: float
    g2: float
    g12: float
    b12: float
    e1: float
    e2: float
    pm1: float
    pm2: float
    h1: float
    h2: float
    d1: float
    d2: float
    fault_g1: float | None = None
    fault_g2: float | None = None


def electrical_power(x: np.ndarray, p: TwoMachineParameters,
                     fault: bool = False) -> Tuple[float, float]:
    """计算两机电磁功率 Pe1、Pe2。

    使用方法：
        传入状态 ``x``（首元素为相对角）和参数，返回 (Pe1, Pe2);``fault=True`` 时
        用故障电导给出故障期功率。
    """

    delta = float(x[0])
    if fault:
        g1 = p.g1 if p.fault_g1 is None else p.fault_g1
        g2 = p.g2 if p.fault_g2 is None else p.fault_g2
        return p.e1 * p.e1 * g1, p.e2 * p.e2 * g2
    pe1 = p.e1 ** 2 * (p.g12 + p.g1) + p.e1 * p.e2 * p.g12 * np.cos(delta) + p.e1 * p.e2 * p.b12 * np.sin(delta)
    pe2 = p.e2 ** 2 * (p.g12 + p.g2) + p.e1 * p.e2 * p.g12 * np.cos(delta) - p.e1 * p.e2 * p.b12 * np.sin(delta)
    return float(pe1), float(pe2)


def f_2m(x: np.ndarray, p: TwoMachineParameters) -> np.ndarray:
    """完整三状态两机摆动方程右端。

    使用方法：
        传入状态 ``x=[delta12, omega12, omega_sum]`` 和参数，返回三维导数向量。
    """

    delta12, omega12, omega_sum = np.asarray(x, dtype=float)
    pe1, pe2 = electrical_power(x, p)
    ddelta = omega12
    domega = p.pm1 / p.h1 - p.pm2 / p.h2 - pe1 / p.h1 + pe2 / p.h2 - omega12 * (p.d1 / p.h1 + p.d2 / p.h2) / 2.0 - omega_sum * (p.d1 / p.h1 - p.d2 / p.h2) / 2.0
    dsum = p.pm1 / p.h1 + p.pm2 / p.h2 - pe1 / p.h1 - pe2 / p.h2 - omega12 * (p.d1 / p.h1 - p.d2 / p.h2) / 2.0 - omega_sum * (p.d1 / p.h1 + p.d2 / p.h2) / 2.0
    return np.array([ddelta, domega, dsum])


def f_2m_fault(x: np.ndarray, p: TwoMachineParameters) -> np.ndarray:
    """故障期间两机摆动方程右端。

    使用方法：
        与 ``f_2m`` 同接口，但电磁功率取故障值，用于故障持续段积分。
    """

    delta12, omega12, omega_sum = np.asarray(x, dtype=float)
    pe1, pe2 = electrical_power(x, p, fault=True)
    domega = p.pm1 / p.h1 - p.pm2 / p.h2 - pe1 / p.h1 + pe2 / p.h2 - omega12 * (p.d1 / p.h1 + p.d2 / p.h2) / 2.0 - omega_sum * (p.d1 / p.h1 - p.d2 / p.h2) / 2.0
    dsum = p.pm1 / p.h1 + p.pm2 / p.h2 - pe1 / p.h1 - pe2 / p.h2 - omega12 * (p.d1 / p.h1 - p.d2 / p.h2) / 2.0 - omega_sum * (p.d1 / p.h1 + p.d2 / p.h2) / 2.0
    return np.array([omega12, domega, dsum])


def f_2m_reduce(x: np.ndarray, p: TwoMachineParameters) -> np.ndarray:
    """去掉共同速度后的二阶相对运动模型右端。

    使用方法：
        传入状态 ``x=[delta12, omega12]`` 和参数，返回二维导数向量，用于相平面分析。
    """

    delta12, omega12 = np.asarray(x, dtype=float)
    pe1, pe2 = electrical_power(np.array([delta12, omega12, 0.0]), p)
    return np.array([
        omega12,
        p.pm1 / p.h1 - p.pm2 / p.h2 - pe1 / p.h1 + pe2 / p.h2 - omega12 * (p.d1 / p.h1 + p.d2 / p.h2) / 2.0,
    ])


def equilibria(p: TwoMachineParameters, guesses: np.ndarray | None = None,
               tolerance: float = 1e-8) -> list[dict]:
    """扫描初值求平衡点并按稳定性分类。

    使用方法：
        传入参数，可选初值网格 ``guesses``，返回平衡点列表；每项含解 x、雅可比、
        特征值和非负特征值个数 ``unstable_dimension``。
    步骤：
        对每个初值用牛顿法求根，去重后线性化求特征值并统计非负实部个数。
    """

    if guesses is None:
        guesses = np.linspace(-2 * np.pi, 2 * np.pi, 17)[:, None]
        guesses = np.c_[guesses, np.zeros((guesses.shape[0], 2))]
    result = []
    for guess in guesses:
        sol, ok, _, _ = newton_solve(lambda z: f_2m(z, p), guess,
                                     tol=tolerance, max_iter=100)
        if not ok:
            continue
        if any(np.linalg.norm(sol - old["x"]) < 10 * tolerance for old in result):
            continue
        jac = numerical_jacobian(lambda z: f_2m(z, p), sol)
        eig = np.linalg.eigvals(jac)
        result.append({"x": sol, "jacobian": jac, "eigenvalues": eig,
                       "unstable_dimension": int(np.sum(np.real(eig) >= 0))})
    return result


def simulate_two_machine(p: TwoMachineParameters, x0: np.ndarray,
                         tlength: float = 10.0, tunit: float = 1e-3,
                         fault_until: float | None = None) -> Tuple[np.ndarray, np.ndarray]:
    """积分两机模型轨迹。

    使用方法：
        传入参数、初值 ``x0``、总时长、步长和故障结束时刻 ``fault_until``，返回
        (时间, 状态);``fault_until`` 之前用故障方程，之后用正常方程。
    """

    n = max(2, int(round(tlength / tunit)))
    time = np.arange(n) * tunit
    states = np.zeros((n, 3), dtype=float)
    states[0] = x0
    for k in range(n - 1):
        fun = f_2m_fault if fault_until is not None and time[k] < fault_until else f_2m
        # 标准 RK4，单独写出四个斜率以便读者追踪。
        k1 = fun(states[k], p)
        k2 = fun(states[k] + 0.5 * tunit * k1, p)
        k3 = fun(states[k] + 0.5 * tunit * k2, p)
        k4 = fun(states[k] + tunit * k3, p)
        states[k + 1] = states[k] + tunit * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
    return time, states
