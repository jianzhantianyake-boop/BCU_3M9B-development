# -*- coding: utf-8 -*-
"""SPM(结构保持模型)能量法: 网络代数解 + 势能 + 发电机功率(对齐 MATLAB Cal_MM_CCT_SPM)。

背景:
    reduced 能量法把负荷母线消去(Kron 约简), SPM 保留网络节点为代数节点。MATLAB 平台的 SPM 用
    **恒阻抗负荷**(并入 Yfull_mod, Sload 的 P/Q=0), 网络母线满足 P=Q=0 的代数约束, 势能含网络电压项。

本模块(逐块对 MATLAB 交叉验证, 均已通过):
    - solve_spm_network: 给定发电机内角, 解网络母线角/电压(恒阻抗, P=Q=0)。对 SEP 网络态 9e-15。
    - spm_generator_power: 发电机经完整网络(Yfull_mod)注入的电磁功率。SEP 处 COI 功率失配 1.9e-11。
    - spm_potential_energy: 5 项 SPM 势能(Ep1 磁势能 / Ep2 电纳 / Ep3 电导 / Ep4 网损 / Ep5 负荷)。
      用 MATLAB 的 CUEP 网络态得 E_crit=3.3757, 对 MATLAB 误差 0。
    - spm_fault_energy_cct: 沿故障轨迹逐点算 Ek+Ep, 首次越过临界能量即 LEA CCT。给定 E_crit=3.3757
      时精确复现 MATLAB 的 CCT=0.2053。

⚠️ 未自足闭环(明确的下一里程碑): 独立求 SPM CUEP 网络态需选对物理分支——SPM 网络方程在 CUEP
发电机角处有 11+ 个解, MATLAB 靠 MGP(Fun_Cal_MGP_SPM/AEiteration_SPM)沿物理轨迹连续跟踪播种。
这是与 v1 find_mgp 同源的数值分支难题, 尚未移植稳健的 SPM 版 controlling-UEP。因此 energy CCT 目前
需外部传入 E_critical(如 3.3757)。物理约束: 正确的 E_crit 须 < 故障能量峰值(本例 5.568), 否则能量
法给不出有限 CCT。

单位: 角度 rad, 功率/导纳 pu。发电机在前 ngen 节点, 网络在后。依赖 scipy。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


def _yfull_mod(state) -> np.ndarray:
    return np.asarray(state.metadata["yfull_mod"], dtype=complex)


# ------------------------- 网络代数解(恒阻抗, P=Q=0) -------------------------

def spm_network_residual(x: np.ndarray, delta_gen: np.ndarray, yfull: np.ndarray,
                         epu: np.ndarray, sload_pq: Optional[np.ndarray] = None) -> np.ndarray:
    """SPM 网络母线功率残差 [P_net; Q_net](恒阻抗负荷 -> 注入为 0; 复现 Fun_AEfslove_SPM)。

    使用方法:
        x=[delta_net(nnet); V_net(nnet)], delta_gen 发电机内角, yfull 重排完整导纳(Yfull_mod),
        epu 发电机内电势。sload_pq(可选)=每网络母线 [P,Q] 附加注入(恒阻抗时为 0)。返回 2*nnet 残差。
    """

    ngen = int(delta_gen.size)
    nbus = int(yfull.shape[0])
    nnet = nbus - ngen
    G = yfull.real
    B = yfull.imag
    dn = x[:nnet]
    vn = x[nnet:]
    P = np.zeros(nnet)
    Q = np.zeros(nnet)
    for i in range(nnet):
        for j in range(ngen):
            dd = dn[i] - delta_gen[j]
            P[i] += vn[i] * epu[j] * (B[i + ngen, j] * np.sin(dd) + G[i + ngen, j] * np.cos(dd))
            Q[i] += vn[i] * epu[j] * (-B[i + ngen, j] * np.cos(dd) + G[i + ngen, j] * np.sin(dd))
        for l in range(nnet):
            dd = dn[i] - dn[l]
            P[i] += vn[i] * vn[l] * (B[i + ngen, l + ngen] * np.sin(dd) + G[i + ngen, l + ngen] * np.cos(dd))
            Q[i] += vn[i] * vn[l] * (-B[i + ngen, l + ngen] * np.cos(dd) + G[i + ngen, l + ngen] * np.sin(dd))
        if sload_pq is not None:
            P[i] += sload_pq[i, 0]
            Q[i] += sload_pq[i, 1]
    return np.r_[P, Q]


def solve_spm_network(delta_gen: np.ndarray, yfull: np.ndarray, epu: np.ndarray,
                      guess: Optional[np.ndarray] = None, sload_pq=None,
                      tol: float = 1e-12, voltage_min: float = 1e-4) -> Tuple[np.ndarray, bool, float]:
    """解 SPM 网络代数方程。

    返回 ``(x=[delta_net; V_net], converged, residual_norm)``。除残差外，
    ``converged`` 还要求所有网络电压大于 ``voltage_min``，以排除方程的
    零电压数学根；连续 warm-start 仍由调用者负责提供物理分支初值。
    """

    from scipy.optimize import root

    ngen = int(delta_gen.size)
    nnet = int(yfull.shape[0]) - ngen
    if guess is None:
        guess = np.r_[np.zeros(nnet), np.ones(nnet)]
    sol = root(lambda x: spm_network_residual(x, delta_gen, yfull, epu, sload_pq),
               np.asarray(guess, dtype=float), method="hybr", tol=tol)
    r = float(np.linalg.norm(spm_network_residual(sol.x, delta_gen, yfull, epu, sload_pq)))
    # The algebraic equations admit mathematically valid zero-voltage roots
    # (especially from cold/random starts).  They are not admissible network
    # states and must never be allowed to masquerade as a physical branch.
    nnet = int(yfull.shape[0]) - ngen
    voltages = np.asarray(sol.x[nnet:], dtype=float)
    physical_voltage = bool(np.all(np.isfinite(voltages)) and
                            np.all(voltages > float(voltage_min)))
    return sol.x, bool(sol.success and r < 1e-6 and physical_voltage), r


# ------------------------- 发电机电磁功率(经 Yfull_mod) -------------------------

def spm_generator_power(delta_gen: np.ndarray, delta_net: np.ndarray, v_net: np.ndarray,
                        yfull: np.ndarray, epu: np.ndarray) -> np.ndarray:
    """发电机经完整网络注入的电磁功率 Pe(用 Yfull_mod, 发电机 V=E)。

    Pe_i = Σ_k E_i V_k (G_ik cos(δ_i-θ_k) + B_ik sin(δ_i-θ_k)), k 遍历全母线(发电机 V=E,θ=δ)。
    """

    ngen = int(delta_gen.size)
    theta = np.r_[delta_gen, delta_net]
    volt = np.r_[epu, v_net]
    G = yfull.real
    B = yfull.imag
    pe = np.zeros(ngen)
    for i in range(ngen):
        dd = theta[i] - theta
        pe[i] = epu[i] * np.sum(volt * (G[i, :] * np.cos(dd) + B[i, :] * np.sin(dd)))
    return pe


# ------------------------- SPM 势能(5 项) -------------------------

def spm_potential_energy(preset, postfault, yfull, sep_gen, sep_net_theta, sep_net_v,
                         end_gen, end_net_theta, end_net_v, sload_full=None) -> np.ndarray:
    """SPM 势能 5 项 [Ep1..Ep5](复现 Fun_Cal_PotentialEnergy_SPM, path_energy_cal=0 时 Ep4=0)。

    使用方法:
        传入 SEP 与末端(通常 CUEP)的发电机角/网络角/网络电压, 返回 [Ep1,Ep2,Ep3,Ep4,Ep5];
        临界能量 = 之和。sload_full(可选)=每网络母线 [P,Q](恒阻抗时全 0, Ep5=0)。
    """

    ngen = int(np.asarray(preset.m).size)
    nbus = int(yfull.shape[0])
    nnet = nbus - ngen
    G = yfull.real
    B = yfull.imag
    Pm = np.asarray(preset.pmpu, dtype=float)
    E = np.asarray(epu_of(preset), dtype=float)

    ep1 = float(np.sum([-(Pm[i] - E[i] ** 2 * G[i, i]) * (end_gen[i] - sep_gen[i]) for i in range(ngen)]))
    ep2 = 0.0
    for i in range(nnet):
        for j in range(ngen):
            ep2 += -(end_net_v[i] * E[j] * B[i + ngen, j] * np.cos(end_net_theta[i] - end_gen[j])
                     - sep_net_v[i] * E[j] * B[i + ngen, j] * np.cos(sep_net_theta[i] - sep_gen[j]))
        ep2 += -(end_net_v[i] ** 2 / 2 * B[i + ngen, i + ngen] - sep_net_v[i] ** 2 / 2 * B[i + ngen, i + ngen])
    for i in range(nnet - 1):
        for l in range(i + 1, nnet):
            ep2 += -(end_net_v[i] * end_net_v[l] * B[i + ngen, l + ngen] * np.cos(end_net_theta[i] - end_net_theta[l])
                     - sep_net_v[i] * sep_net_v[l] * B[i + ngen, l + ngen] * np.cos(sep_net_theta[i] - sep_net_theta[l]))
    ep3 = float(np.sum([G[i + ngen, i + ngen] / 3 * (end_net_theta[i] - sep_net_theta[i])
                        * (end_net_v[i] ** 2 + end_net_v[i] * sep_net_v[i] + sep_net_v[i] ** 2) for i in range(nnet)]))
    ep4 = 0.0  # path_energy_cal=0: Ray 近似为 0
    ep5 = 0.0
    if sload_full is not None:
        for i in range(nnet):
            ep5 += sload_full[i, 0] * (end_net_theta[i] - sep_net_theta[i]) \
                + sload_full[i, 1] * (np.log(end_net_v[i]) - np.log(sep_net_v[i]))
    return np.array([ep1, ep2, ep3, ep4, float(ep5)])


def epu_of(preset) -> np.ndarray:
    return np.asarray(preset.epu, dtype=float)


# ------------------------- SPM 能量法 CCT -------------------------

def spm_fault_energy_cct(static, e_critical: float, tfault: float = 0.6,
                         tunit: float = 1e-4) -> Tuple[float, bool]:
    """SPM 能量法 CCT: 沿故障轨迹逐点算 Ek+Ep(postfault 势能), 首次越 e_critical 即 CCT。

    使用方法:
        传入 v1 StaticResult 与临界能量 e_critical(SPM 势能意义, 如与 MATLAB 交叉验证得 3.3757),
        返回 (CCT[s], 是否找到)。给定 e_critical=3.3757 时对 MATLAB 的 0.2053 精确复现。
    机理(复现 Fun_Cal_CCT_Energy_SPM):
        (1) 发电机角轨迹用 reduced 故障积分(Kron 约简对恒阻抗负荷精确, 故 SPM 发电机角=reduced);
        (2) 每步在 **postfault 网络** 上解网络代数(warm-start, 分支连续)得网络角/电压;
        (3) Ek=0.5 Σ m ωc², Ep=spm_potential_energy(SEP->当前), 首次 Ek+Ep>e_critical 的时刻即 CCT。
    诚实说明: e_critical 需外部提供(自足求 CUEP 网络态的分支选择尚未闭环, 见模块头部)。
    """

    from bcu_3m9b.dynamics import integrate_reduced

    preset, base, post, fault = static.preset, static.basevalue, static.postfault, static.fault
    ypost = _yfull_mod(post)
    epu = np.asarray(preset.epu, dtype=float)
    m = np.asarray(preset.m, dtype=float)
    ngen = int(m.size)
    nnet = int(ypost.shape[0]) - ngen

    # SEP 参考网络态.
    xsep, ok, _ = solve_spm_network(np.asarray(post.sep_delta, dtype=float), ypost, epu)
    sep_gen = np.asarray(post.sep_delta, dtype=float)
    sep_nt, sep_nv = xsep[:nnet], xsep[nnet:]

    # 故障网络母线号 -> postfault 网络索引(能量在 postfault 全维度上算, 故障母线补 0).
    faultbus = int(preset.fault_line[preset.fault_position])
    post_net = [b for b in np.asarray(post.metadata["transform"]).astype(int) if b > ngen]
    fpos = post_net.index(faultbus)
    yfault = _yfull_mod(fault)
    nnet_f = int(yfault.shape[0]) - ngen

    d0 = np.asarray(static.prefault.sep_delta, dtype=float)
    w0 = np.full(ngen, static.prefault.sep_omegapu * base.omega_b)
    traj = integrate_reduced(tfault, tunit, fault, preset, base, d0, w0)

    guess = xsep.copy()
    e_prev = None
    for k in range(traj.time.size):
        thc = traj.thetac[k]
        x, _, _ = solve_spm_network(thc, ypost, epu, guess=guess)
        guess = x
        nt, nv = x[:nnet], x[nnet:]
        ek = 0.5 * float(np.sum(m * traj.omegac[k] ** 2))
        ep = float(np.sum(spm_potential_energy(preset, post, ypost, sep_gen, sep_nt, sep_nv, thc, nt, nv)))
        esum = ek + ep
        if e_prev is not None and e_prev < e_critical < esum:
            return float((k - 1) * tunit), True
        e_prev = esum
    return float(traj.time[-1]), False
