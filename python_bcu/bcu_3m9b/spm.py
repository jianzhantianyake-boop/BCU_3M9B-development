"""结构保持模型（SPM）教学实现使用说明。

使用方法：
    给定发电机角度求负荷母线代数解用 ``solve_algebraic``；算发电机导数用
    ``spm_generator_rhs``；跑一段 SPM 显式 DAE 近似实验用 ``simulate_spm``。

设计说明：
    原 MATLAB 版本用带 Mass 矩阵的 ``ode15s`` 解 DAE。这里用同一功率平衡方程，
    但每个显式时间步先用纯 NumPy 牛顿法求负荷母线电压，再推进发电机摆动方程。
    这是可实际运行的实验实现，不伪装成 MATLAB 求解器的逐点复现。
"""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np

from .equilibrium import electrical_power
from .numerics import newton_solve
from .types import BaseValue, NetworkState, Preset


def _load_power(preset: Preset, bus_numbers: np.ndarray) -> np.ndarray:
    """按重排后的普通母线编号取 P/Q 负荷标幺值。

    使用方法：
        传入参数和普通母线编号序列，返回形状 (n, 2) 的 [P, Q]；无负荷数据时返回零。
    """

    result = np.zeros((bus_numbers.size, 2), dtype=float)
    if preset.s_load is None:
        return result
    for k, bus in enumerate(bus_numbers):
        row = np.flatnonzero(preset.s_load[:, 0].astype(int) == int(bus))
        if row.size:
            result[k] = preset.s_load[row[0], 1:3]
    return result


def algebraic_residual(z: np.ndarray, delta_gen: np.ndarray,
                       yorg_ordered: np.ndarray, load_pq: np.ndarray,
                       ngen: int) -> np.ndarray:
    """计算负荷母线的有功/无功平衡残差。

    使用方法：
        传入负荷代数量 ``z=[角度; 幅值]``、发电机角度、重排导纳、负荷 P/Q 和发电机
        台数，返回残差向量，配合牛顿法求负荷节点电压。
    """

    nload = yorg_ordered.shape[0] - ngen
    delta_load = z[:nload]
    voltage_load = z[nload:]
    voltage = np.r_[np.exp(1j * delta_gen), voltage_load * np.exp(1j * delta_load)]
    injection = voltage * np.conj(yorg_ordered @ voltage)
    # 负荷为正消耗，故网络注入 + 负荷 = 0。
    return np.r_[injection.real[ngen:] + load_pq[:, 0],
                 injection.imag[ngen:] + load_pq[:, 1]]


def solve_algebraic(delta_gen: np.ndarray, state: NetworkState,
                    preset: Preset, guess: np.ndarray | None = None,
                    tol: float = 1e-8) -> Tuple[np.ndarray, bool, float]:
    """给定发电机内部角度，求负荷节点角度与电压。

    使用方法：
        传入发电机角度、含重排信息的网络工况和参数，返回 (代数解 z, 是否收敛, 残差)；
        ``guess`` 可传上一步解作为热启动加速收敛。
    说明：
        求解顺序: scipy.optimize.root(hybr, 信赖域, 对冷启动更稳) 依次尝试[热启动 guess ->
        平启动 -> 若干扰动初值], 任一达标即返回; 全部不理想时回退 v1 自写阻尼牛顿, 取残差最小
        者。这样根治原 v1 隐患(纯 numpy 牛顿在故障后代数解冷启动偶发首步不收敛)。需要严格 DAE
        级稳健化(刚性隐式积分 + 连续法)见 v2 的 spm_dae.py。
    """

    from scipy.optimize import root

    yorg = np.asarray(state.metadata.get("yorg_mod", state.yfull), dtype=complex)
    transform = np.asarray(state.metadata.get("transform"), dtype=int)
    ngen = preset.ngen
    load_pq = _load_power(preset, transform[ngen:])
    nload = yorg.shape[0] - ngen

    def resid(z):
        return algebraic_residual(z, delta_gen, yorg, load_pq, ngen)

    default = np.r_[np.zeros(nload), np.ones(nload)]
    threshold = max(1e3 * tol, 1e-8)
    candidates = []
    if guess is not None:
        candidates.append(np.asarray(guess, dtype=float))
    candidates.append(default)
    # 冷启动兜底: 负荷角小扰动 + 电压略降的几个物理合理初值.
    rng = np.array([0.1, -0.1, 0.2])
    for a in rng:
        candidates.append(np.r_[np.full(nload, a), np.full(nload, 0.95)])

    best, best_r = default, np.inf
    for g in candidates:
        sol = root(resid, g, method="hybr", tol=tol)
        r = float(np.linalg.norm(resid(sol.x)))
        if r < best_r:
            best, best_r = sol.x, r
        if sol.success and r < threshold:
            return sol.x, True, r
    # scipy 都不达标时回退 v1 阻尼牛顿(仍取全局残差最小者返回).
    sol2, _, _, r2 = newton_solve(resid, candidates[0], tol=tol, max_iter=100, jac_step=1e-6)
    if r2 < best_r:
        best, best_r = sol2, r2
    return best, bool(best_r < threshold), best_r


def spm_generator_rhs(delta_gen: np.ndarray, omega: np.ndarray,
                      state: NetworkState, preset: Preset,
                      basevalue: BaseValue, algebraic: np.ndarray) -> np.ndarray:
    """用负荷母线解算发电机角度和速度导数。

    使用方法：
        传入发电机角度、角速度、网络工况、参数、基值和当前代数解 ``algebraic``，
        返回 [dtheta; domega]；电磁功率由完整网络注入的实部得到。
    """

    yorg = np.asarray(state.metadata.get("yorg_mod", state.yfull), dtype=complex)
    transform = np.asarray(state.metadata.get("transform"), dtype=int)
    nload = yorg.shape[0] - preset.ngen
    voltage = np.r_[np.exp(1j * delta_gen), algebraic[nload:] * np.exp(1j * algebraic[:nload])]
    injection = voltage * np.conj(yorg @ voltage)
    pe = injection.real[:preset.ngen]
    coi = np.dot(omega, preset.m) / np.sum(preset.m)
    return np.r_[omega - coi,
                 (preset.pmpu - pe - preset.d * (omega - basevalue.omega_b)) / preset.m]


def simulate_spm(tlength: float, tunit: float, state: NetworkState,
                 preset: Preset, basevalue: BaseValue,
                 delta0: np.ndarray, omega0: np.ndarray,
                 guess: np.ndarray | None = None) -> Dict[str, np.ndarray]:
    """跑一段结构保持模型的显式 DAE 近似。

    使用方法：
        传入总时长、步长、网络工况、参数、基值和初值，返回字典 {time, delta,
        omega, algebraic, delta_coi, omega_coi}。
    步骤：
        每步先解负荷代数方程（热启动），再用 RK4 推进发电机摆动方程；代数方程不
        收敛会抛出异常并报告所在步与残差。
    """

    n = preset.ngen
    steps = max(2, int(round(tlength / tunit)))
    time = np.arange(steps) * tunit
    delta = np.zeros((steps, n))
    omega = np.zeros((steps, n))
    algebraic = []
    delta[0], omega[0] = delta0, omega0
    # guess: 第一步负荷代数解的热启动初值(缺省 None 用内部 [0..., 1...] 默认);
    # 后续步自动用上一步解热启动.
    for k in range(steps):
        z, ok, residual = solve_algebraic(delta[k], state, preset, guess)
        if not ok:
            raise RuntimeError(f"SPM load algebraic equations did not converge; step={k}, residual={residual:g}")
        guess = z
        algebraic.append(z.copy())
        if k == steps - 1:
            break
        x = np.r_[delta[k], omega[k]]
        fun = lambda xx: spm_generator_rhs(xx[:n], xx[n:], state, preset, basevalue, z)
        k1 = fun(x)
        k2 = fun(x + 0.5 * tunit * k1)
        k3 = fun(x + 0.5 * tunit * k2)
        k4 = fun(x + tunit * k3)
        xn = x + tunit * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        delta[k + 1], omega[k + 1] = xn[:n], xn[n:]
    return {"time": time, "delta": delta, "omega": omega,
            "algebraic": np.asarray(algebraic),
            "delta_coi": delta - (delta @ preset.m / np.sum(preset.m))[:, None],
            "omega_coi": omega - (omega @ preset.m / np.sum(preset.m))[:, None]}
