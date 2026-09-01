# -*- coding: utf-8 -*-
"""P2.2 / P2.4: 单机级参考动态模型(可验证的基础, 尚未接入 BCU/能量函数流水线).

内容:
    P2.2 one-axis(磁链衰减)同步机 SMIB 模型: 状态 [δ, ω, E'q], 含 q 轴暂态 EMF 动态与励磁 Efd.
    P2.4 GFM(构网型)下垂逆变器 SMIB 模型: P-ω 下垂, 状态 [δ](可扩 Q-V).

定位(诚实说明):
    这些是**单机无穷大母线级**的正确参考实现, 在极限/稳态下可校验(见 test_models.py), 用作后续
    多机化与能量函数扩展的基础; 目前**未**接入 BCU/CCT/稳定域流水线(那是更大的工作).
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

OMEGA_S = 2.0 * np.pi * 60.0


# ============================ P2.2: one-axis 同步机 ============================

@dataclass
class OneAxisSMIB:
    """one-axis(磁链衰减)同步机接无穷大母线参数.

    参数: H(s), D, Xd, Xd1(=Xd'), Xq, Xe(外接电抗), Vinf(无穷大母线电压), Pm, Efd, Tdo1(=T'do,s).
    round-rotor 取 Xq=Xd1.
    """

    H: float
    D: float
    Xd: float
    Xd1: float
    Xq: float
    Xe: float
    Vinf: float
    Pm: float
    Efd: float
    Tdo1: float
    omega_s: float = OMEGA_S


def one_axis_currents(delta, Eq1, p: OneAxisSMIB):
    """由 (δ, E'q) 求 d/q 轴电流 (Id, Iq)(经无穷大母线闭式)."""

    Id = (Eq1 - p.Vinf * np.cos(delta)) / (p.Xd1 + p.Xe)
    Iq = (p.Vinf * np.sin(delta)) / (p.Xq + p.Xe)
    return Id, Iq


def one_axis_power(delta, Eq1, p: OneAxisSMIB) -> float:
    """电磁功率 Pe = E'q Iq + (Xq - Xd') Id Iq(凸极); Xq=Xd' 时退化为 E'q Iq."""

    Id, Iq = one_axis_currents(delta, Eq1, p)
    return float(Eq1 * Iq + (p.Xq - p.Xd1) * Id * Iq)


def one_axis_rhs(x, p: OneAxisSMIB):
    """one-axis 模型右端: 状态 x=[δ, ω(rad/s slip 绝对角速度), E'q], 返回导数.

    dδ/dt = ω - ωs; dω/dt = (ωs/2H)(Pm - Pe - D(ω-ωs)); dE'q/dt = (Efd - E'q - (Xd-Xd')Id)/T'do.
    """

    delta, omega, Eq1 = float(x[0]), float(x[1]), float(x[2])
    Id, Iq = one_axis_currents(delta, Eq1, p)
    Pe = Eq1 * Iq + (p.Xq - p.Xd1) * Id * Iq
    ddelta = omega - p.omega_s
    domega = (p.omega_s / (2.0 * p.H)) * (p.Pm - Pe - p.D * (omega - p.omega_s))
    dEq1 = (p.Efd - Eq1 - (p.Xd - p.Xd1) * Id) / p.Tdo1
    return np.array([ddelta, domega, dEq1])


# ============================ P2.4: GFM 下垂逆变器 ============================

@dataclass
class GFMDroopSMIB:
    """构网型(GFM)下垂逆变器接无穷大母线参数.

    参数: E(逆变器电压幅值 pu), Vinf, X(并网电抗), Pset(有功设定 pu), mp(P-ω 下垂系数),
    omega_s.
    功率-功角: P(δ)=E*Vinf*sin(δ)/X; 一阶 P-ω 下垂: dδ/dt = ωs*mp*(Pset - P).
    """

    E: float
    Vinf: float
    X: float
    Pset: float
    mp: float
    omega_s: float = OMEGA_S


def gfm_power(delta, p: GFMDroopSMIB) -> float:
    """GFM 输出有功 P(δ)=E*Vinf*sin(δ)/X."""

    return float(p.E * p.Vinf * np.sin(delta) / p.X)


def gfm_rhs(x, p: GFMDroopSMIB):
    """GFM 下垂右端: 状态 x=[δ], dδ/dt = ωs*mp*(Pset - P(δ))."""

    delta = float(x[0])
    return np.array([p.omega_s * p.mp * (p.Pset - gfm_power(delta, p))])


def gfm_equilibrium(p: GFMDroopSMIB) -> float:
    """GFM 稳态功角 δss = asin(Pset*X/(E*Vinf))(存在则返回, 否则 ValueError)."""

    s = p.Pset * p.X / (p.E * p.Vinf)
    if abs(s) > 1.0:
        raise ValueError("no GFM equilibrium (Pset too large for the line)")
    return float(np.arcsin(s))
