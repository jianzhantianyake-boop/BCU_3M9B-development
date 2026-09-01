# -*- coding: utf-8 -*-
"""SPM(结构保持模型)能量法: 网络代数解 + 势能 + 发电机功率(对齐 MATLAB Cal_MM_CCT_SPM)。

背景:
    reduced 能量法把负荷母线消去(Kron 约简), SPM 保留网络节点为代数节点。MATLAB 平台的 SPM 用
    **恒阻抗负荷**(并入 Yfull_mod, Sload 的 P/Q=0), 网络母线满足 P=Q=0 的代数约束, 势能含网络电压项。

本模块的底层公式曾做过分块 MATLAB 对照；当前报告只把这些结果作为带限制的历史证据：
    - solve_spm_network: 给定发电机内角, 解网络母线角/电压(恒阻抗, P=Q=0)。对 SEP 网络态 9e-15。
    - spm_generator_power: 发电机经完整网络(Yfull_mod)注入的电磁功率。SEP 处 COI 功率失配 1.9e-11。
    - spm_potential_energy: 5 项 SPM 势能(Ep1 磁势能 / Ep2 电纳 / Ep3 电导 / Ep4 网损 / Ep5 负荷)。
      MATLAB 曾输出 E_crit=3.3757，但其 CUEP 网络角与 raw fsolve 角存在混合坐标，当前不作为
      物理一致性证据。
    - spm_fault_energy_cct: 沿实际 fault-network DAE 轨迹逐点算 Ek+Ep，首次越过临界能量即 LEA CCT；
      旧的 E_crit=3.3757 -> 0.2053 仅是历史代理流程，不代表当前严格 DAE 结果。

⚠️ 未自足闭环(明确的下一里程碑): 独立求 SPM CUEP 网络态需选对物理分支——SPM 网络方程在 CUEP
发电机角处有 11+ 个解, MATLAB 靠 MGP(Fun_Cal_MGP_SPM/AEiteration_SPM)沿物理轨迹连续跟踪播种。
这是与 v1 find_mgp 同源的数值分支难题, 尚未移植稳健的 SPM 版 controlling-UEP。因此 energy CCT 目前
需外部传入 E_critical(如 3.3757)。物理约束: 正确的 E_crit 须 < 故障能量峰值；当前严格
fault-network DAE 轨迹在默认案例上的峰值约为 4.82202，旧的 5.568/6.46886 仅是历史代理值；
若临界能量不低于该峰值，能量法就给不出有限 CCT。

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

def _expand_fault_network_state(z_fault: np.ndarray, fault, postfault,
                                ngen: int) -> tuple[np.ndarray, np.ndarray]:
    """Map a fault-network algebraic state into the postfault bus ordering.

    The MATLAB SPM trajectory keeps a zero placeholder for the removed fault
    bus.  Python's physical fault solve omits that bus, so the mapping is done
    by recorded bus numbers rather than a hard-coded column index.
    """
    z_fault = np.asarray(z_fault, dtype=float).reshape(-1)
    n_fault = int(z_fault.size // 2)
    theta_fault = z_fault[:n_fault]
    voltage_fault = z_fault[n_fault:]
    fault_transform = np.asarray(fault.metadata["transform"], dtype=int)[ngen:]
    post_transform = np.asarray(postfault.metadata["transform"], dtype=int)[ngen:]
    if fault_transform.size != n_fault:
        raise ValueError("fault algebraic state width does not match fault transform")
    theta_post = np.zeros(post_transform.size, dtype=float)
    voltage_post = np.zeros(post_transform.size, dtype=float)
    for k, bus in enumerate(fault_transform):
        matches = np.flatnonzero(post_transform == int(bus))
        if matches.size != 1:
            raise ValueError(f"fault bus {int(bus)} is missing or duplicated in postfault transform")
        j = int(matches[0])
        theta_post[j] = theta_fault[k]
        voltage_post[j] = voltage_fault[k]
    return theta_post, voltage_post


def spm_fault_energy_series(static, *, tfault: float = 0.6,
                            tunit: float = 1e-4, method: str = "Radau",
                            max_points: int | None = None) -> tuple[np.ndarray, np.ndarray, bool]:
    """Compute total SPM energy on the actual fault-network DAE trajectory.

    The differential trajectory and algebraic states both come from
    :func:`bcu_v2.spm_dae.simulate_spm_dae`.  The postfault network and SEP are
    used only as the potential-energy reference, matching the MATLAB energy
    functional.  A false return flag means that no complete finite series is
    available; callers must not replace it with a reduced-model result.
    """
    if tfault <= 0 or tunit <= 0:
        raise ValueError("tfault and tunit must be positive")
    from .spm_dae import simulate_spm_dae

    preset, base = static.preset, static.basevalue
    fault, postfault = static.fault, static.postfault
    ypost = _yfull_mod(postfault)
    epu = np.asarray(preset.epu, dtype=float)
    ngen = int(preset.ngen)
    nnet = int(ypost.shape[0]) - ngen
    sep_gen = np.asarray(postfault.sep_delta, dtype=float)
    sep_state, sep_ok, _ = solve_spm_network(sep_gen, ypost, epu)
    if not sep_ok:
        return np.array([], dtype=float), np.array([], dtype=float), False
    sep_nt, sep_nv = sep_state[:nnet], sep_state[nnet:]
    delta0 = np.asarray(static.prefault.sep_delta, dtype=float)
    omega0 = np.full(ngen, float(static.prefault.sep_omegapu) * base.omega_b)
    trajectory = simulate_spm_dae(tfault, tunit, fault, preset, base,
                                  delta0, omega0, method=method)
    if not trajectory.get("success", False):
        return np.asarray(trajectory.get("time", []), dtype=float), \
            np.full(len(trajectory.get("time", [])), np.nan), False

    times = np.asarray(trajectory["time"], dtype=float)
    count = times.size
    if max_points is not None:
        if max_points <= 0:
            raise ValueError("max_points must be positive when provided")
        indices = np.linspace(0, max(count - 1, 0), min(max_points, count)).astype(int)
        indices = np.unique(indices)
    else:
        indices = np.arange(count, dtype=int)
    energies = np.full(indices.size, np.nan, dtype=float)
    m = np.asarray(preset.m, dtype=float)
    for out_index, k in enumerate(indices):
        dg = np.asarray(trajectory["delta"][k], dtype=float)
        omega = np.asarray(trajectory["omega"][k], dtype=float)
        z_fault = np.asarray(trajectory["algebraic"][k], dtype=float)
        if not (np.all(np.isfinite(dg)) and np.all(np.isfinite(omega))
                and np.all(np.isfinite(z_fault))):
            continue
        theta_net, voltage_net = _expand_fault_network_state(
            z_fault, fault, postfault, ngen
        )
        if not (np.all(np.isfinite(theta_net)) and np.all(np.isfinite(voltage_net))):
            continue
        omega_coi = omega - np.dot(m, omega) / np.sum(m)
        ep = spm_potential_energy(
            preset, postfault, ypost, sep_gen, sep_nt, sep_nv,
            dg, theta_net, voltage_net,
        )
        energies[out_index] = float(np.sum(ep) + 0.5 * np.sum(m * omega_coi ** 2))
    return times[indices], energies, bool(energies.size and np.all(np.isfinite(energies)))

def spm_fault_energy_cct(static, e_critical: float, tfault: float = 0.6,
                         tunit: float = 1e-4) -> Tuple[float, bool]:
    """SPM 能量法 CCT: 沿 fault-network DAE 逐点算 Ek+Ep, 首次越界即 CCT。

    使用方法:
        传入 v1 StaticResult 与临界能量 e_critical(SPM 势能意义), 返回 (CCT[s], 是否找到)。
        历史 MATLAB 代理值不再被视为当前严格 DAE 的验收基准。
    机理(保留 MATLAB 势能泛函):
        (1) fault-only SPM DAE 同时推进发电机状态并连续校正故障网络代数状态;
        (2) 将删去的故障母线按记录的母线号映射为零占位，再用 postfault 网络计算势能;
        (3) Ek=0.5 Σ m ωc², Ep=spm_potential_energy(SEP->当前), 首次 Ek+Ep>e_critical 的时刻即 CCT。
    诚实说明: e_critical 仍需调用者提供；自足 CUEP 由 ``spm_cuep`` 负责计算。
    """
    if not np.isfinite(e_critical):
        return float("nan"), False
    times, energies, valid = spm_fault_energy_series(
        static, tfault=tfault, tunit=tunit, method="Radau", max_points=None,
    )
    if not valid or times.size == 0:
        return float("nan"), False
    previous = float(energies[0])
    if previous >= e_critical:
        return float(times[0]), True
    for k in range(1, times.size):
        current = float(energies[k])
        if previous < e_critical <= current:
            return float(times[k - 1]), True
        previous = current
    return float(times[-1]), False
