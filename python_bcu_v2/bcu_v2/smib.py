# -*- coding: utf-8 -*-
"""单机无穷大母线(SMIB) + 等面积准则闭式 CCT(验证金标准).

使用方法:
    用 SMIB 描述"发电机-无穷大母线"三阶段(故障前/故障中/故障后)功角特性; 等面积准则
    给出临界切除角 δcc 的闭式解, 对 Pmax_fault=0(母线金属性短路)还给出临界切除时间 tcc
    的闭式解. 这两个闭式解与 cct.py 里事件驱动的数值 CCT 相互印证, 是不依赖任何实现的
    正确性金标准.

约定:
    状态 x=[δ(rad), ν(rad/s)], ν=dδ/dt. 摆动方程(标幺, 惯性常数 H):
        dδ/dt = ν
        dν/dt = (Pm - Pmax_stage*sin(δ) - D*ν) / M,   M = 2H/ωs
    等面积闭式解要求无阻尼 D=0.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass
class SMIB:
    """SMIB 参数容器.

    参数:
        H: 惯性常数(s). D: 阻尼系数(标幺, 对 ν). Pm: 机械功率(pu).
        Pmax_pre/fault/post: 三阶段的功率-功角幅值 Pmax=E*V/X(pu).
        omega_s: 基准角速度(rad/s), 默认 60 Hz.
    """

    H: float
    Pm: float
    Pmax_pre: float
    Pmax_fault: float
    Pmax_post: float
    D: float = 0.0
    omega_s: float = 2.0 * np.pi * 60.0

    @property
    def M(self) -> float:
        # 用法: 惯性 M=2H/ωs.
        return 2.0 * self.H / self.omega_s


def delta0(s: SMIB) -> float:
    """故障前稳定平衡角 δ0 = asin(Pm/Pmax_pre)."""

    return float(np.arcsin(s.Pm / s.Pmax_pre))


def delta_sep_post(s: SMIB) -> float:
    """故障后稳定平衡角 δs = asin(Pm/Pmax_post)."""

    return float(np.arcsin(s.Pm / s.Pmax_post))


def delta_uep_post(s: SMIB) -> float:
    """故障后不稳定平衡角(上界) δu = π - asin(Pm/Pmax_post)."""

    return float(np.pi - np.arcsin(s.Pm / s.Pmax_post))


def critical_angle(s: SMIB) -> float:
    """等面积准则闭式临界切除角 δcc.

    使用方法: 传入 SMIB, 返回临界切除角(rad). 公式:
        cos δcc = [Pm(δu-δ0) - Pmax_f cosδ0 + Pmax_p cosδu] / (Pmax_p - Pmax_f)
    若右端绝对值 >1(该工况无临界角/恒稳或恒不稳)抛出 ValueError.
    """

    d0 = delta0(s)
    du = delta_uep_post(s)
    denom = s.Pmax_post - s.Pmax_fault
    if abs(denom) < 1e-12:
        raise ValueError("Pmax_post == Pmax_fault, equal-area formula degenerate")
    c = (s.Pm * (du - d0) - s.Pmax_fault * np.cos(d0) + s.Pmax_post * np.cos(du)) / denom
    if abs(c) > 1.0:
        raise ValueError(f"no critical angle (cos={c:.4g}); system always stable or always unstable")
    return float(np.arccos(c))


def critical_time_analytic(s: SMIB) -> float:
    """临界切除时间闭式解(仅当 Pmax_fault=0 且 D=0).

    使用方法: 母线金属性短路(故障期功率为 0)且无阻尼时, 故障期匀加速:
        δ(t) = δ0 + 0.5*(Pm/M)*t^2  =>  tcc = sqrt(2M(δcc-δ0)/Pm).
    其他情形无闭式 tcc(需数值积分故障段), 本函数抛出 ValueError.
    """

    if abs(s.Pmax_fault) > 1e-12 or abs(s.D) > 1e-12:
        raise ValueError("analytic tcc only valid for Pmax_fault=0 and D=0")
    dcc = critical_angle(s)
    return float(np.sqrt(2.0 * s.M * (dcc - delta0(s)) / s.Pm))


def swing_rhs(x: np.ndarray, s: SMIB, stage: str) -> np.ndarray:
    """SMIB 摆动方程右端.

    使用方法: 传入状态 x=[δ, ν], SMIB, 阶段 stage in {'pre','fault','post'}, 返回 [dδ, dν].
    """

    pmax = {"pre": s.Pmax_pre, "fault": s.Pmax_fault, "post": s.Pmax_post}[stage]
    delta, nu = float(x[0]), float(x[1])
    return np.array([nu, (s.Pm - pmax * np.sin(delta) - s.D * nu) / s.M])


def default_smib(bolted: bool = True) -> SMIB:
    """返回一个自洽的示例 SMIB.

    使用方法: bolted=True 时 Pmax_fault=0(便于同时验证闭式 tcc); False 时含转移功率.
    """

    return SMIB(H=3.5, Pm=1.0, Pmax_pre=2.0,
                Pmax_fault=0.0 if bolted else 0.5,
                Pmax_post=1.5, D=0.0)
