# -*- coding: utf-8 -*-
"""P2.1: 通用系统装配层(任意 MATPOWER case + 每机动态参数), 并配置 39 母线动态案例.

使用方法:
    build_preset(case, H, Xd1, ...) 把"每台机的惯性常数 H、暂态电抗 Xd'"等按 case 的发电机
    顺序装成 v1 的 Preset(m=2H/ωs), 从而让 v1 的 build_static_result 能跑**任意** case, 不再
    局限 9 母线. case39_dynamic() 返回配好动态参数的 IEEE 39 母线(10 机)案例.

单位约定(与 v1 一致):
    m = 2H/ωs(pu-s^2), d 为阻尼系数, pmpu=Pg/Sbase, xd1=Xd'(pu), epu=发电机内电势(经典模型
    flag_xd=0 时取端电压设定值 Vg).

重要提示:
    39 母线的 H、Xd' 为**示例值**(接近标准 New England 数据), 定量 CCT 前请用你的权威参考
    (如 Pai 1989 附录 / Athay 1979)核对替换. 装配机理与 SEP/潮流不依赖这些数值的精确性.
"""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np

OMEGA_S = 2.0 * np.pi * 60.0  # 基准角速度 rad/s


def build_preset(case, H, Xd1, epu=None, damping=0.1,
                 fault_line=None, fault_position: int = 0, flag_xd: int = 0):
    """把每机动态参数按 case 的发电机顺序装成 v1 Preset.

    参数:
        case: CaseData(其 gen 行顺序决定各数组顺序).
        H: 每机惯性常数(s), 数组, 与 case.gen 行对齐. m 由 2H/ωs 得到.
        Xd1: 每机暂态电抗 Xd'(pu), 数组.
        epu: 每机内电势(pu); 缺省取 case.gen 的 Vg 列(经典模型 flag_xd=0).
        damping: 标量(则 d=damping*m)或每机数组.
        fault_line: 故障线路 [i, j]; 缺省用 v1 默认 [9,6].
        fault_position: 故障母线在 fault_line 中的索引.
    返回:
        v1 的 Preset(可直接送入 build_static_result).
    """

    from bcu_3m9b.types import Preset

    gen = np.asarray(case.gen, dtype=float)
    H = np.asarray(H, dtype=float)
    Xd1 = np.asarray(Xd1, dtype=float)
    m = 2.0 * H / OMEGA_S
    d = (damping * m) if np.isscalar(damping) else np.asarray(damping, dtype=float)
    pmpu = gen[:, 1] / case.base_mva          # Pg / Sbase
    epu = np.asarray(epu, dtype=float) if epu is not None else gen[:, 5].copy()  # Vg
    kw = {}
    if fault_line is not None:
        kw["fault_line"] = np.asarray(fault_line, dtype=int)
    return Preset(m=m, d=d, pmpu=pmpu, xd1=Xd1, epu=epu,
                  flag_xd=flag_xd, fault_position=fault_position, **kw)


# ----------------- IEEE 39 母线(10 机)示例动态数据 -----------------
# 键为发电机所在母线号(30..39). H(s)、Xd'(pu), 100 MVA 基值. **示例值, 请核对替换**.
NE39_H = {30: 42.0, 31: 30.3, 32: 35.8, 33: 28.6, 34: 26.0,
          35: 34.8, 36: 26.4, 37: 24.3, 38: 34.5, 39: 500.0}
NE39_XD1 = {30: 0.031, 31: 0.0697, 32: 0.0531, 33: 0.0436, 34: 0.132,
            35: 0.05, 36: 0.049, 37: 0.057, 38: 0.057, 39: 0.006}


def case39_dynamic(fault_line=(16, 17), fault_position: int = 0):
    """返回配好动态参数的 IEEE 39 母线(10 机)案例.

    使用方法: 返回 (case, preset); 可送入 build_static_result(case, preset) 做 39 母线静态初始化.
    注意: H/Xd' 为示例值(见模块顶部提示).
    """

    from bcu_3m9b.cases import case39_modified

    case = case39_modified()
    buses = case.gen[:, 0].astype(int)  # 发电机母线顺序
    H = np.array([NE39_H[b] for b in buses], dtype=float)
    Xd1 = np.array([NE39_XD1[b] for b in buses], dtype=float)
    preset = build_preset(case, H, Xd1, fault_line=fault_line, fault_position=fault_position)
    return case, preset


def build_static_dynamic(case, preset, solve_tol: float = 1e-9):
    """对任意(case, preset)做静态初始化(潮流->约简->SEP), 复用 v1.build_static_result.

    使用方法: 返回 v1 的 StaticResult; 是"把 9 母线专用流程推广到任意系统"的入口.
    """

    from bcu_3m9b.bcu import build_static_result
    return build_static_result(case=case, preset=preset, solve_tol=solve_tol)
