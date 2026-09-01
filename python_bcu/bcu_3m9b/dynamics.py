"""约简模型暂态轨迹、退出点和时域 CCT 使用说明。

使用方法：
    需要右端函数用 ``reduced_rhs``；定步长积分一段轨迹用 ``integrate_reduced``；
    找退出点用 ``find_exitpoint``；判轨迹是否稳定用 ``trajectory_stable``；分段
    仿真故障+切除用 ``simulate_fault_clear``；用固定网格估时域 CCT 用
    ``time_domain_cct``。
"""

from __future__ import annotations

from typing import Callable, Dict, Optional, Tuple

import numpy as np

from .equilibrium import electrical_power
from .types import BaseValue, NetworkState, Preset, Trajectory


def reduced_rhs(x: np.ndarray, yred: np.ndarray, preset: Preset,
                basevalue: BaseValue) -> np.ndarray:
    """计算网络约简摆动方程的右端项。

    使用方法：
        传入打包状态 ``x = [theta; omega]``（绝对功角与绝对角速度）、约简导纳
        ``yred`` 和参数 ``preset``、``basevalue``，返回 ``[dtheta; domega]``。
    参数：
        x：长度 2n 的状态向量，前 n 个功角（rad），后 n 个角速度（pu）。
        yred：n×n 约简导纳矩阵。
        preset：机组惯量 m、阻尼 d、机械功率 pmpu、内电势 epu 等。
        basevalue：系统基值，omega_b 为基准角速度（rad/s）。
    返回：
        长度 2n 的导数向量 [dtheta; domega]。
    步骤：
        (1) 拆分状态为功角与角速度；(2) 计算电磁功率；(3) 从角速度减去 COI 漂移；
        (4) 由 机械功率 - 电磁功率 - 阻尼 得到加速度。
    """

    n = preset.ngen
    theta = np.asarray(x[:n], dtype=float)
    omega = np.asarray(x[n:2 * n], dtype=float)
    pe = electrical_power(theta, yred, preset.epu)
    coi = np.dot(omega, preset.m) / np.sum(preset.m)
    dtheta = omega - coi
    domega = (preset.pmpu - pe - preset.d * (omega - basevalue.omega_b)) / preset.m
    return np.r_[dtheta, domega]


def integrate_reduced(tlength: float, tunit: float, state: NetworkState,
                      preset: Preset, basevalue: BaseValue,
                      delta0: np.ndarray, omega0: np.ndarray,
                      semi_rk4: bool = True) -> Trajectory:
    """定步长积分一段约简模型轨迹。

    使用方法：
        传入总时长 ``tlength``、步长 ``tunit``、网络工况和初值 ``delta0``、
        ``omega0``，返回 ``Trajectory``。``semi_rk4=True`` 保留原 MATLAB
        Fun_TrajIter_SRF 的计算顺序：电磁功率在一个时间步内固定为步首值，只有
        阻尼项使用中间速度；关闭时改用标准四阶 RK4，便于对方法本身做教学对比。
    参数：
        tlength：积分总时长（s）。
        tunit：固定步长（s）。
        state：提供 ``yred`` 的网络工况。
        delta0/omega0：初始功角与角速度。
    返回：
        含时间、绝对量、COI 相对量和逐步电磁功率的 ``Trajectory``。
    """

    n = preset.ngen
    cycle = max(2, int(round(tlength / tunit)))
    theta = np.zeros((cycle, n), dtype=float)
    omega = np.zeros((cycle, n), dtype=float)
    thetac = np.zeros_like(theta)
    omegac = np.zeros_like(omega)
    omegacoi = np.zeros(cycle, dtype=float)
    pe_all = np.zeros_like(theta)
    theta[0] = np.asarray(delta0, dtype=float)
    omega[0] = np.asarray(omega0, dtype=float)
    for k in range(cycle):
        # 每步先记录电磁功率和 COI 相对量，供能量法和退出点判据复用。
        pe = electrical_power(theta[k], state.yred, preset.epu)
        pe_all[k] = pe
        coi_theta = np.dot(theta[k], preset.m) / np.sum(preset.m)
        coi_omega = np.dot(omega[k], preset.m) / np.sum(preset.m)
        thetac[k] = theta[k] - coi_theta
        omegac[k] = omega[k] - coi_omega
        omegacoi[k] = coi_omega
        if k == cycle - 1:
            break
        if semi_rk4:
            # 复刻原 MATLAB 文件中 S_ 的四列公式：功率固定为步首值，只有阻尼用中间速度。
            accel = np.empty((n, 4), dtype=float)
            stage_fraction = (0.0, 0.5, 0.5, 1.0)
            for j in range(4):
                omega_mid = omega[k] + stage_fraction[j] * tunit * (
                    accel[:, j - 1] if j else 0.0
                )
                accel[:, j] = (preset.pmpu - pe - preset.d * (omega_mid - basevalue.omega_b)) / preset.m
            speed = np.column_stack([
                omega[k], omega[k] + 0.5 * tunit * accel[:, 0],
                omega[k] + 0.5 * tunit * accel[:, 1],
                omega[k] + tunit * accel[:, 2],
            ])
            theta[k + 1] = theta[k] + (speed @ np.array([1.0, 2.0, 2.0, 1.0])) * tunit / 6.0
            omega[k + 1] = omega[k] + (accel @ np.array([1.0, 2.0, 2.0, 1.0])) * tunit / 6.0
        else:
            # 标准 RK4：四个斜率均重算电磁功率，逐个写出便于对照教学。
            x = np.r_[theta[k], omega[k]]
            k1 = reduced_rhs(x, state.yred, preset, basevalue)
            k2 = reduced_rhs(x + 0.5 * tunit * k1, state.yred, preset, basevalue)
            k3 = reduced_rhs(x + 0.5 * tunit * k2, state.yred, preset, basevalue)
            k4 = reduced_rhs(x + tunit * k3, state.yred, preset, basevalue)
            xn = x + tunit * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
            theta[k + 1], omega[k + 1] = xn[:n], xn[n:]
    time = np.arange(cycle, dtype=float) * tunit
    return Trajectory(time, theta, omega, thetac, omegac, omegacoi, pe_all, tunit)


def find_exitpoint(traj: Trajectory, postfault: NetworkState,
                   preset: Preset, eps: float = 1e-9) -> int:
    """在轨迹上寻找退出点。

    使用方法：
        传入故障轨迹和故障后工况，返回退出点的时间索引；判据为功率失配与相对速度点积由
        "明显为负"变"明显为正"的首个过零点。
    说明：
        排除初始 ``ωc≈0`` 造成的伪过零(原 v1 隐患: 默认故障下退出点退化为 index=0); 用 ``eps``
        阈值确保过零点两侧点积有明确符号, 并跳过第 0 个点。
    """

    values = np.zeros(traj.time.size)
    for k in range(traj.time.size):
        pe = electrical_power(traj.theta[k], postfault.yred, preset.epu)
        values[k] = np.dot(preset.pmpu - pe, traj.omegac[k])
    mask = (values[:-1] < -eps) & (values[1:] > eps)
    if mask.size:
        mask[0] = False  # 排除初始伪过零
    crossing = np.flatnonzero(mask)
    return int(crossing[0]) if crossing.size else traj.time.size - 1


def trajectory_stable(traj: Trajectory, postfault: NetworkState,
                      preset: Preset, tail_tolerance: float = 0.1,
                      criterion: str = "bounded", bound: float = 1.5 * np.pi) -> bool:
    """判断轨迹是否稳定。

    使用方法：
        传入轨迹和故障后工况，返回 True/False。总是先查失步硬判据(任意两机角差≥2π 即失稳);
        再按 ``criterion`` 选稳定判据。
    参数：
        criterion: "bounded"(默认, 有界性判据) 或 "return_to_sep"(原 v1 "末端回到 SEP" 判据)。
        bound: 有界性判据的允许偏离(相对故障后 SEP), 默认 1.5π。
        tail_tolerance: return_to_sep 判据的末端容差。
    说明：
        默认改用**有界性判据**(全程有限且离故障后 SEP 不超过 bound): 原 v1 的"末端回到 SEP"
        判据在**轻阻尼**系统下失效(阻尼小, 末端仍在 SEP 附近振荡而非收敛, 被误判失稳 -> REA≈0)。
        需完全复刻原行为时传 criterion="return_to_sep"。
    """

    if np.max(np.abs(traj.theta[:, :, None] - traj.theta[:, None, :])) >= 2.0 * np.pi:
        return False
    if criterion == "bounded":
        if postfault.sep_delta is None:
            return bool(np.all(np.isfinite(traj.thetac)))
        return bool(np.all(np.isfinite(traj.thetac))
                    and np.max(np.abs(traj.thetac - postfault.sep_delta)) < bound)
    if criterion == "return_to_sep":
        if postfault.sep_delta is None:
            raise ValueError("return_to_sep criterion requires the post-fault SEP")
        return bool(np.linalg.norm(traj.thetac[-1] - postfault.sep_delta) <= tail_tolerance)
    raise ValueError(f"unknown criterion {criterion!r}")


def simulate_fault_clear(clear_time: float, fault_duration: float,
                         postfault_duration: float, tunit: float,
                         fault: NetworkState, postfault: NetworkState,
                         preset: Preset, basevalue: BaseValue,
                         delta0: np.ndarray, omega0: np.ndarray) -> Dict[str, object]:
    """做一次“故障—切除”的分段仿真。

    使用方法：
        传入切除时刻 ``clear_time``、各段时长、故障与故障后工况和初值，返回字典
        {fault, postfault, stable}；故障段积分到切除时刻，用末端状态作为故障后段
        初值继续积分并判稳。
    """

    fault_traj = integrate_reduced(clear_time, tunit, fault, preset, basevalue, delta0, omega0)
    delta_clear = fault_traj.theta[-1]
    omega_clear = fault_traj.omega[-1]
    post_traj = integrate_reduced(postfault_duration, tunit, postfault, preset, basevalue,
                                  delta_clear, omega_clear)
    return {
        "fault": fault_traj,
        "postfault": post_traj,
        "stable": trajectory_stable(post_traj, postfault, preset),
    }


def time_domain_cct(fault_duration: float, postfault_duration: float,
                    tunit: float, fault: NetworkState, postfault: NetworkState,
                    preset: Preset, basevalue: BaseValue,
                    delta0: np.ndarray, omega0: np.ndarray,
                    samples: int = 21) -> Tuple[float, Dict[str, object]]:
    """用固定网格搜索时域 CCT（REA 的可运行平台入口）。

    使用方法：
        传入各段时长、工况、初值和网格点数 ``samples``，返回 (最大稳定切除时间,
        细节字典)；在 [0, fault_duration] 上等距取切除时刻，逐个判稳，遇到首个失
        稳即停并返回上一个稳定时刻。
    """

    clear_times = np.linspace(0.0, fault_duration, samples)
    last = None
    stable_time = 0.0
    for clear_time in clear_times:
        current = simulate_fault_clear(clear_time, fault_duration, postfault_duration,
                                        tunit, fault, postfault, preset, basevalue,
                                        delta0, omega0)
        last = current
        if current["stable"]:
            stable_time = float(clear_time)
        else:
            break
    return stable_time, {"last": last, "tested_clear_times": clear_times}
