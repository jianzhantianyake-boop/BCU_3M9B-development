"""MATLAB 同名函数门面使用说明。

使用方法：
    平台内部用更清晰的模块化接口；本文件提供与 ``Fun_*`` 同名的薄封装，便于把
    MATLAB 文件逐个对照阅读。所有函数只用显式传入的参数，不读全局变量。
"""

from __future__ import annotations

import numpy as np

from .dynamics import find_exitpoint, integrate_reduced, trajectory_stable
from .energy import energy_cct, find_mgp, potential_energy
from .equilibrium import sep_check, sep_residual, solve_sep
from .network import kron_reduce, rxb_to_yfull
from .types import BaseValue, CCTResult, NetworkState, Preset, Trajectory


def Fun_RXB2Yfull(RXB, pfdata):
    """Fun_RXB2Yfull 门面：由 RXB 线路参数拼完整导纳。

    使用方法：传入 RXB 矩阵和 pfdata，返回 ``Yfull``。
    """

    return rxb_to_yfull(RXB, pfdata.nbus)


def Fun_Yfull2Yred(Y_full, pfdata, faultflag=0):
    """Fun_Yfull2Yred 门面：对完整导纳做 Kron 约简。

    使用方法：传入 ``Y_full``、pfdata 和故障标志；``faultflag`` 为含故障母线的序
    列时会先删除该母线，返回分块与约简结果。
    """

    active = np.arange(1, Y_full.shape[0] + 1, dtype=int)
    removed = None
    if isinstance(faultflag, (list, tuple, np.ndarray)) and len(faultflag) > 1:
        removed = int(faultflag[1])
        active = np.asarray([b for b in range(1, pfdata.nbus + 1) if b != removed])
    return kron_reduce(Y_full, pfdata.gen_no, active, removed)


def Fun_SEPfslove(delta_omega, preset, state, basevalue):
    """Fun_SEPfslove 门面：返回 SEP 残差函数值。

    使用方法：传入打包量和参数，返回 n 维残差，供求根使用。
    """

    return sep_residual(delta_omega, preset, state, basevalue)


def Fun_SEPcheck(state, preset, delta_sep, omega_sep):
    """Fun_SEPcheck 门面：核对 SEP 残差和偏离标志。

    使用方法：传入工况、参数、待检角度与速度，返回 (残差, 是否偏离)。
    """

    return sep_check(state, preset, delta_sep, omega_sep)


def Fun_SEPiteration(yred, pmpu, epu, m, d, delta0, omega0, omegab,
                     n_itermax=10000, tolerr=1e-8):
    """Fun_SEPiteration 门面：用统一牛顿求解器替代手写迭代。

    使用方法：传入约简导纳和机组参数及初值，返回 (COI 角度, pu 速度, 收敛标志0/1,
    迭代次数)。
    """

    preset = Preset(np.asarray(m), np.asarray(d), np.asarray(pmpu),
                    np.zeros_like(np.asarray(m)), np.asarray(epu))
    state = NetworkState("compat", np.asarray(yred), np.asarray(yred),
                         np.empty((0, 0)), np.empty((0, 0)), np.empty((0, 0)), np.empty((0, 0)))
    base = BaseValue(omegab)
    delta, omegapu, _, ok, iterations = solve_sep(preset, state, base,
                                                   np.asarray(delta0), omega0,
                                                   tolerr, n_itermax)
    return delta, omegapu, int(ok), iterations


def Fun_TrajIter_SRF(tlength, tunit, yred, preset, delta0, omega0, omegab):
    """Fun_TrajIter_SRF 门面：按 MATLAB 返回顺序给出轨迹分量。

    使用方法：传入时长、步长、约简导纳、参数和初值，返回
    (theta, omega, thetac, omegacoi, pe, 步数)。
    """

    state = NetworkState("compat", np.asarray(yred), np.asarray(yred),
                         np.empty((0, 0)), np.empty((0, 0)), np.empty((0, 0)), np.empty((0, 0)))
    traj = integrate_reduced(tlength, tunit, state, preset, BaseValue(omegab), delta0, omega0)
    return traj.theta, traj.omega, traj.thetac, traj.omegacoi, traj.pe, traj.time.size


def Fun_Cal_Exitpoint(tlength, tunit, yredfault, yredpost, preset, delta0, omega0, omegab):
    """Fun_Cal_Exitpoint 门面：积分故障轨迹并返回退出点索引。

    使用方法：传入时长、步长、故障与故障后约简导纳、参数和初值，返回
    (theta, omega, thetac, omegac, 退出点索引)。
    """

    fault = NetworkState("fault", yredfault, yredfault, *(np.empty((0, 0)) for _ in range(4)))
    post = NetworkState("postfault", yredpost, yredpost, *(np.empty((0, 0)) for _ in range(4)))
    traj = integrate_reduced(tlength, tunit, fault, preset, BaseValue(omegab), delta0, omega0)
    index = find_exitpoint(traj, post, preset)
    return traj.theta, traj.omega, traj.thetac, traj.omegac, index


def Fun_Cal_PotentialEnergy(preset, postfault, thetac_start, thetac_end):
    """Fun_Cal_PotentialEnergy 门面：返回三项势能。

    使用方法：传入参数、故障后工况和路径起止角，返回 (Ep1, Ep2, Ep3)。
    """

    return tuple(potential_energy(preset, postfault, thetac_start, thetac_end))


def Fun_Cal_MGP(thetac_escape, postfault, preset):
    """Fun_Cal_MGP 门面：非交互式 MGP 搜索。

    使用方法：传入退出点角、故障后工况和参数，返回
    (MGP 角, 1, 是否找到, 残差范数序列, 命中残差)。
    """

    result = find_mgp(thetac_escape, postfault, preset)
    return result["theta_mgp"], 1, int(result["found"]), result["norm"], float(result["norm"][result["index"]])


def Fun_Cal_CCT_Energy(critical_energy, fault_traj, postfault, preset):
    """Fun_Cal_CCT_Energy 门面：能量法求 CCT。

    使用方法：传入临界能量、故障轨迹、故障后工况和参数，返回
    (CCT, thetac, omegac, 是否越界)。
    """

    result = energy_cct(critical_energy, fault_traj, preset, postfault)
    return result.cct, result.exit_state.get("thetac"), result.exit_state.get("omegac"), result.flag_cct


def Fun_TrajIter_StableCheck_SRF(tlength, tunit, postfault, preset, delta0, omega0, omegab):
    """Fun_TrajIter_StableCheck_SRF 门面：积分并给稳定性标志。

    使用方法：传入时长、步长、故障后工况、参数和初值，返回
    (theta, omega, thetac, omegacoi, pe, 步数, 不稳定标志)。
    """

    traj = integrate_reduced(tlength, tunit, postfault, preset, BaseValue(omegab), delta0, omega0)
    return traj.theta, traj.omega, traj.thetac, traj.omegacoi, traj.pe, traj.time.size, int(not trajectory_stable(traj, postfault, preset))
