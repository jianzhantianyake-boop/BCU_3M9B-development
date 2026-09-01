# -*- coding: utf-8 -*-
"""P1.2: 事件驱动 + 二分搜索的精确临界切除时间(CCT).

使用方法:
    用 scipy.solve_ivp 的自适应步长 + 事件检测精确定位"失步时刻", 再对切除时刻做二分,
    得到比固定网格精确 1~2 个数量级的 CCT. 提供 SMIB 版(可与等面积闭式解印证)与 3 机
    网络约简版(与 v1 的网格 REA / 能量法 LEA 对比).

依赖: scipy.
"""

from __future__ import annotations

from typing import Callable, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from . import smib as _smib


def bisect_cct(is_stable: Callable[[float], bool], hi: float = 0.5,
               lo: float = 0.0, tol: float = 1e-5, max_expand: int = 25) -> Tuple[float, bool]:
    """对"切除时刻->是否稳定"做二分, 返回 (CCT, 是否找到边界).

    使用方法:
        传入判稳函数 is_stable(tc); 假设 tc=lo 稳定. 先把 hi 上探到不稳定, 再二分到区间
        小于 tol. 若始终稳定则返回当前 hi 与 False.
    """

    if not is_stable(lo):
        return lo, False
    h = hi
    for _ in range(max_expand):
        if not is_stable(h):
            break
        h *= 1.5
    else:
        return h, False  # 探测范围内始终稳定
    hi = h
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if is_stable(mid):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), True


# ------------------------------- SMIB -------------------------------

def precise_cct_smib(s: _smib.SMIB, t_post: float = 5.0, tol: float = 1e-5,
                     rtol: float = 1e-10, atol: float = 1e-12) -> Tuple[float, bool]:
    """SMIB 的事件驱动精确 CCT.

    使用方法:
        对候选切除时刻 tc: 故障段积分到 tc, 故障后段积分并以"δ 越过故障后 UEP"为失步事件;
        无事件即稳定. 二分得 CCT. 可与 smib.critical_time_analytic / critical_angle 对比.
    """

    d0 = _smib.delta0(s)
    du = _smib.delta_uep_post(s)

    def unstable_event(t, z):
        return z[0] - du
    unstable_event.terminal = True
    unstable_event.direction = 1.0

    def is_stable(tc: float) -> bool:
        if tc <= 0:
            return True
        s1 = solve_ivp(lambda t, z: _smib.swing_rhs(z, s, "fault"),
                       [0.0, tc], [d0, 0.0], rtol=rtol, atol=atol)
        xc = s1.y[:, -1]
        if xc[0] >= du:
            return False  # 切除时 δ 已越过故障后 UEP -> 必失步(事件检测不到起点越界)
        s2 = solve_ivp(lambda t, z: _smib.swing_rhs(z, s, "post"),
                       [0.0, t_post], xc, events=unstable_event, rtol=rtol, atol=atol)
        return len(s2.t_events[0]) == 0

    return bisect_cct(is_stable, hi=0.3, tol=tol)


# --------------------------- 3 机网络约简 ---------------------------

def precise_cct_reduced(static, t_post: float = 3.0, tol: float = 1e-4,
                        sep_limit: float = 2.0 * np.pi,
                        rtol: float = 1e-8, atol: float = 1e-10) -> Tuple[float, bool]:
    """网络约简模型的事件驱动精确 CCT.

    使用方法:
        以预故障 SEP 为初值; 故障段积分到 tc, 故障后段以"任意两机角差超过 sep_limit"为
        失步事件, 无事件且有界即稳定; 二分得 CCT. 与 v1 的网格 REA、能量法 LEA 对照.
    """

    from bcu_3m9b.dynamics import reduced_rhs

    preset, base = static.preset, static.basevalue
    n = preset.ngen
    d0 = np.asarray(static.prefault.sep_delta, dtype=float)
    w0 = np.full(n, static.prefault.sep_omegapu * base.omega_b)
    x0 = np.r_[d0, w0]

    def rhs_of(state_net):
        return lambda t, x: reduced_rhs(x, state_net.yred, preset, base)

    def sep_event(t, x):
        th = x[:n]
        return (np.max(th) - np.min(th)) - sep_limit
    sep_event.terminal = True
    sep_event.direction = 1.0

    def is_stable(tc: float) -> bool:
        if tc <= 0:
            return True
        s1 = solve_ivp(rhs_of(static.fault), [0.0, tc], x0, rtol=rtol, atol=atol)
        xc = s1.y[:, -1]
        s2 = solve_ivp(rhs_of(static.postfault), [0.0, t_post], xc,
                       events=sep_event, rtol=rtol, atol=atol)
        if len(s2.t_events[0]) > 0:
            return False
        return bool(np.all(np.isfinite(s2.y[:, -1])))

    return bisect_cct(is_stable, hi=0.3, tol=tol)
