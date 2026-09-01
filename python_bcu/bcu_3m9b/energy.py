"""势能、动能、MGP 和能量法 CCT 使用说明。

使用方法：
    算三项势能用 ``potential_energy``；算动能用 ``kinetic_energy``；对整条轨迹算
    Ep/Ek/总能量用 ``trajectory_energy``；用能量法找 CCT 用 ``energy_cct``；追踪
    MGP 单轨迹用 ``mgp_single_trajectory``，取 MGP 候选点用 ``find_mgp``。
"""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np

from .equilibrium import electrical_power
from .types import CCTResult, NetworkState, Preset, Trajectory


def potential_energy(preset: Preset, postfault: NetworkState,
                    theta_start: np.ndarray, theta_end: np.ndarray) -> np.ndarray:
    """计算三项势能分解（复现 Fun_Cal_PotentialEnergy）。

    使用方法：
        传入参数、故障后工况、起点角 ``theta_start`` 和终点角 ``theta_end``，
        返回长度 3 的 [Ep1, Ep2, Ep3]；总势能为三者之和。
    参数：
        theta_start：路径起点功角（COI 坐标）。
        theta_end：路径终点功角（COI 坐标）。
    步骤：
        Ep1 为磁势能相对 SEP 的位移项；Ep2 为电纳耦合的余弦项；Ep3 为电导耦合的
        路径相关项，按 ``preset.path_energy_cal`` 选择忽略/解析近似/多段梯形积分。
    """

    sep = np.asarray(postfault.sep_delta, dtype=float)
    start = np.asarray(theta_start, dtype=float).reshape(-1)
    end = np.asarray(theta_end, dtype=float).reshape(-1)
    g = np.real(postfault.yred)
    b = np.imag(postfault.yred)
    e = preset.epu
    n = preset.ngen
    ep1 = -np.sum((preset.pmpu - e * e * np.diag(g)) * (end - sep))
    ep2 = 0.0
    ep3 = 0.0
    for i in range(n - 1):
        for j in range(i + 1, n):
            dij_end = end[i] - end[j]
            dij_sep = sep[i] - sep[j]
            ep2 += -e[i] * e[j] * b[i, j] * (np.cos(dij_end) - np.cos(dij_sep))
            di = end[i] - start[i]
            dj = end[j] - start[j]
            dij = di - dj
            ratio = (di + dj) / dij if abs(dij) > 1e-7 else di + dj
            ep3 += e[i] * e[j] * g[i, j] * ratio * (np.sin(dij_end) - np.sin(start[i] - start[j]))
    if preset.path_energy_cal == -1:
        # path_energy_cal=-1：忽略电导路径项，只保留磁势能与电纳项。
        ep3 = 0.0
    elif preset.path_energy_cal > 0:
        # path_energy_cal>0：教学版多段梯形积分，对每个有序交叉项逐段累加。
        segments = max(1, int(preset.path_energy_cal))
        ep3 = 0.0
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                di = (end[i] - start[i]) / segments
                dj = (end[j] - start[j]) / segments
                for k in range(segments):
                    d0 = (start[i] + k * di) - (start[j] + k * dj)
                    d1 = (start[i] + (k + 1) * di) - (start[j] + (k + 1) * dj)
                    ep3 += 0.5 * e[i] * e[j] * g[i, j] * di * (np.cos(d0) + np.cos(d1))
    return np.array([float(ep1), float(ep2), float(ep3)])


def kinetic_energy(preset: Preset, omegac: np.ndarray) -> float:
    """计算 COI 动能 0.5*m*omega_c^2。

    使用方法：
        传入参数和 COI 相对速度 ``omegac``，返回标量动能。
    """

    return float(0.5 * np.dot(preset.m, np.asarray(omegac) ** 2))


def trajectory_energy(traj: Trajectory, preset: Preset,
                      postfault: NetworkState) -> Dict[str, np.ndarray]:
    """算一条轨迹上的势能、动能和总能量。

    使用方法：
        传入轨迹、参数和故障后工况，返回字典 {ep_components, ep, ek, total}，
        其中 total = ep + ek，供能量法判据使用。
    """

    ep = np.vstack([
        potential_energy(preset, postfault, postfault.sep_delta, theta)
        for theta in traj.thetac
    ])
    ek = np.array([kinetic_energy(preset, speed) for speed in traj.omegac])
    return {"ep_components": ep, "ep": ep.sum(axis=1), "ek": ek, "total": ep.sum(axis=1) + ek}


def energy_cct(critical_energy: float, traj: Trajectory, preset: Preset,
               postfault: NetworkState) -> CCTResult:
    """用能量法找总能量首次越过临界能量的时刻。

    使用方法：
        传入临界能量、故障轨迹、参数和故障后工况，返回 ``CCTResult``；找到首个
        总能量由小于变大于临界值的离散时刻即为 LEA CCT，未越界则返回末端并置
        ``flag_cct=False``。
    """

    data = trajectory_energy(traj, preset, postfault)
    crossing = np.flatnonzero((data["total"][:-1] < critical_energy)
                              & (data["total"][1:] > critical_energy))
    if crossing.size:
        index = int(crossing[0])
        return CCTResult(float(traj.time[index]), index,
                         {"thetac": traj.thetac[index], "omegac": traj.omegac[index]},
                         True, "LEA", {"critical_energy": critical_energy, "energy": data})
    return CCTResult(float(traj.time[-1]), traj.time.size - 1,
                     {"thetac": traj.thetac[-1], "omegac": traj.omegac[-1]},
                     False, "LEA", {"critical_energy": critical_energy, "energy": data})


def mgp_single_trajectory(theta_start: np.ndarray, yred: np.ndarray,
                          preset: Preset, tunit: float = 1e-4,
                          n_itermax: int = 10, norm_tol: float = 1e-5) -> Tuple[np.ndarray, np.ndarray, int, bool]:
    """沿阻尼边界追踪 MGP 单轨迹。

    使用方法：
        传入起点角、约简导纳和参数，返回 (轨迹, 残差范数序列, 命中索引, 是否找到)；
        用于给 CUEP 求解提供初值猜测。
    步骤：
        每步算 COI 功率残差范数并沿 residual/d 方向推进；当范数在低值区出现回升
        （局部极小）时判定命中并返回该点索引。
    """

    n = preset.ngen
    theta = np.zeros((n_itermax, n), dtype=float)
    norms = np.zeros(n_itermax, dtype=float)
    theta[0] = theta_start
    msum = np.sum(preset.m)
    for k in range(n_itermax):
        pe = electrical_power(theta[k], yred, preset.epu)
        pcoi = np.sum(preset.pmpu - pe)
        residual = preset.pmpu - pe - pcoi / msum * preset.m
        norms[k] = np.linalg.norm(residual)
        if k < n_itermax - 1:
            speed_like = residual / preset.d
            next_theta = theta[k] - theta[k, -1]
            next_theta[:-1] += (speed_like[:-1] - speed_like[-1]) * tunit
            next_theta[-1] = 0.0
            theta[k + 1] = next_theta - np.dot(preset.m, next_theta) / msum
        if k >= 2 and norms[k] > norms[k - 1] + norm_tol and norms[k - 1] < 0.1:
            return theta, norms, k, True
    return theta, norms, int(np.argmin(norms)), False


def find_mgp(theta_escape: np.ndarray, postfault: NetworkState,
             preset: Preset, n_itermax: int = 10) -> Dict[str, object]:
    """非交互式 MGP 搜索。

    使用方法：
        传入退出点角 ``theta_escape``、故障后工况和参数，返回字典
        {theta_mgp, trajectory, norm, index, found}；实验平台优先返回可检查的
        候选点，不强行断言严格 MGP。
    """

    theta, norms, index, found = mgp_single_trajectory(
        np.asarray(theta_escape), postfault.yred, preset, n_itermax=n_itermax
    )
    return {
        "theta_mgp": theta[index].copy(),
        "trajectory": theta,
        "norm": norms,
        "index": index,
        "found": found,
    }
