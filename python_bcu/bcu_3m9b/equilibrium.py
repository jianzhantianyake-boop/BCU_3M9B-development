"""SEP/CUEP 残差与平衡点求解使用说明。

使用方法：
    求某台机电磁功率用 ``electrical_power``；把角度转到 COI 坐标用
    ``normalize_coi``；求稳定平衡点 SEP 用 ``solve_sep``（会回填到
    ``NetworkState``）；核对 SEP 用 ``sep_check``；以 MGP 为初值找 CUEP 用
    ``solve_cuep_from_guess``。
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from .numerics import newton_solve
from .types import BaseValue, NetworkState, Preset


def electrical_power(delta: np.ndarray, yred: np.ndarray,
                     epu: np.ndarray) -> np.ndarray:
    """按 BCU 公式计算每台发电机的电磁有功功率。

    使用方法：
        传入功角 ``delta``（rad）、约简导纳 ``yred`` 和内电势 ``epu``，返回各机
        电磁有功向量。
    步骤：
        用角度差矩阵与 G、B 分量按 E_i E_j (G cosΔ + B sinΔ) 求和。
    """

    delta = np.asarray(delta, dtype=float).reshape(-1)
    epu = np.asarray(epu, dtype=float).reshape(-1)
    difference = delta[:, None] - delta[None, :]
    g = np.real(yred)
    b = np.imag(yred)
    return np.sum(
        epu[:, None] * epu[None, :] * (g * np.cos(difference) + b * np.sin(difference)),
        axis=1,
    )


def normalize_coi(delta: np.ndarray, m: np.ndarray) -> np.ndarray:
    """把角度变换到质量加权的 COI 坐标。

    使用方法：
        传入角度向量和惯量 ``m``，返回减去 COI 参考后的角度，使加权和为零。
    """

    delta = np.asarray(delta, dtype=float).reshape(-1)
    m = np.asarray(m, dtype=float).reshape(-1)
    return delta - np.dot(m, delta) / np.sum(m)


def sep_residual(delta_omega: np.ndarray, preset: Preset,
                 state: NetworkState, basevalue: BaseValue) -> np.ndarray:
    """计算 SEP 的 n 维功率残差（对应 Fun_SEPfslove）。

    使用方法：
        传入打包量 ``delta_omega``（前 n-1 个相对角 + 1 个速度）及参数，返回
        n 维残差，配合 ``newton_solve`` 求根。
    """

    z = np.asarray(delta_omega, dtype=float).reshape(-1)
    n = preset.ngen
    delta = np.r_[z[: n - 1], 0.0]
    omega = z[n - 1]
    pe = electrical_power(delta, state.yred, preset.epu)
    pcoi = np.sum(preset.pmpu - pe) - np.sum(preset.d) * omega
    residual = np.empty(n, dtype=float)
    residual[: n - 1] = (
        preset.pmpu[: n - 1] - pe[: n - 1]
        - preset.m[: n - 1] / np.sum(preset.m) * pcoi
        - preset.d[: n - 1] * omega
    )
    residual[n - 1] = np.sum(preset.pmpu - pe) - np.sum(preset.d) * omega
    return residual


def solve_sep(preset: Preset, state: NetworkState, basevalue: BaseValue,
              delta0: Optional[np.ndarray] = None,
              omega0: float = 0.0, tol: float = 1e-9,
              max_iter: int = 100, inplace: bool = True) -> Tuple[np.ndarray, float, np.ndarray, bool, int]:
    """求解稳定平衡点 SEP。

    使用方法：
        传入参数与网络工况和初值，返回 (COI 角度, pu 速度, 功率残差, 收敛标志,
        迭代次数)；``inplace=True``（默认）时把结果回填到 ``state`` 的 sep_* 字段。
    参数：
        delta0：功角初值（缺省取零向量）。
        omega0：速度初值（rad/s）。
        tol：牛顿收敛阈值。
        max_iter：最大迭代步数。
        inplace：是否回填 ``state.sep_*``。静态初始化(build_static_result)需要回填故默认
            True；求 CUEP 等"以异地平衡点为目标"的调用应传 False，避免污染 ``state`` 的 SEP
            字段(原 v1 隐患: 求 CUEP 会改写 postfault.SEP, 跨模式共用 static 出错)。
    步骤：
        先把初值归一到 COI，打包成前 n-1 个相对角加速度求根，收敛后还原 COI
        角度并计算 pu 速度和残差。
    """

    n = preset.ngen
    if delta0 is None:
        delta0 = np.zeros(n, dtype=float)
    delta0 = normalize_coi(delta0, preset.m)
    z0 = np.r_[delta0[: n - 1] - delta0[-1], omega0]
    solution, success, iterations, _ = newton_solve(
        lambda z: sep_residual(z, preset, state, basevalue),
        z0, tol=tol, max_iter=max_iter, jac_step=1e-6,
        name=f"{state.name} SEP",
    )
    delta_raw = np.r_[solution[: n - 1], 0.0]
    delta = normalize_coi(delta_raw, preset.m)
    omega = solution[-1]
    perr = sep_residual(solution, preset, state, basevalue)
    success = bool(success and np.linalg.norm(perr) <= max(10 * tol, 1e-7))
    omegapu = 1.0 + omega / basevalue.omega_b
    if inplace:
        state.sep_delta = delta
        state.sep_omegapu = omegapu
        state.sep_perr = perr
    return delta, omegapu, perr, success, iterations


def sep_check(state: NetworkState, preset: Preset, delta: np.ndarray,
              omega: float) -> Tuple[np.ndarray, bool]:
    """核对给定点是否偏离平衡（对应 Fun_SEPcheck）。

    使用方法：
        传入网络工况、参数、待检角度和速度，返回 (功率残差, 是否明显偏离)；
        判据为残差二范数大于 1e-2。
    """

    pe = electrical_power(delta, state.yred, preset.epu)
    pcoi = np.sum(preset.pmpu) - np.sum(pe) - np.sum(preset.d) * omega
    perr = preset.pmpu - pe - pcoi / np.sum(preset.m) * preset.m - preset.d * omega
    return perr, bool(np.linalg.norm(perr) > 1e-2)


def solve_cuep_from_guess(preset: Preset, state: NetworkState,
                          basevalue: BaseValue, guess: np.ndarray,
                          max_iter: int = 200) -> Tuple[np.ndarray, float, np.ndarray, bool, int]:
    """以 MGP/UEP 猜测为初值搜索控制不稳定平衡点 CUEP。

    使用方法：
        传入参数、故障后工况、基值和猜测点，返回与 ``solve_sep`` 相同的五元组；若牛顿法回落
        到 SEP，调用者可据此改记近似临界点。以 ``inplace=False`` 调用，**不**改写 ``state`` 的
        SEP 字段(避免求 CUEP 污染 postfault.SEP)。
    """

    return solve_sep(preset, state, basevalue, delta0=guess, omega0=0.0,
                     tol=1e-8, max_iter=max_iter, inplace=False)
