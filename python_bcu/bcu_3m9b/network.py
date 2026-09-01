"""网络构造、故障设置和 Kron 约简使用说明。

使用方法：
    先用 ``rxb_to_yfull`` 由线路参数拼出完整导纳，用 ``add_load_admittance``
    加入恒阻抗负荷；正常工况用 ``kron_reduce`` 得到 ``Yred``；故障工况用
    ``remove_fault_bus`` 删除故障母线、切除后用 ``remove_fault_line`` 删线并重拼；
    结构保持模型用 ``reorder_structure_preserved`` 把发电机节点排到前面。

对应关系：
    实现 MATLAB Yfull/Yred 桥梁的等价流程，Kron 消元公式为
    ``Yred = Ynn - Ynr inv(Yrr) Yrn``。
"""

from __future__ import annotations

from typing import Iterable, Optional, Tuple

import numpy as np

from .types import PFData


def rxb_to_yfull(rxb: np.ndarray, nbus: int) -> np.ndarray:
    """由 RXB 线路参数组装完整导纳矩阵（对应 Fun_RXB2Yfull）。

    使用方法：
        传入 [i, j, r, x, b] 线路矩阵和母线数，返回复数 ``Yfull``。
    参数：
        rxb：每行 [起点, 终点, r, x, b]。
        nbus：母线总数。
    步骤：
        逐条线路加串联导纳 1/(r+jx) 到对角和非对角，半充电电容 jb/2 加到两端对角。
    """

    yfull = np.zeros((nbus, nbus), dtype=complex)
    for line in np.asarray(rxb, dtype=float):
        i, j = int(line[0]) - 1, int(line[1]) - 1
        r, x, b = line[2], line[3], line[4]
        den = r * r + x * x
        if den <= 0:
            raise ValueError(f"Line {i + 1}-{j + 1} cannot have both r and x equal to zero")
        y_series = 1.0 / complex(r, x)
        y_shunt = 1j * b / 2.0
        yfull[i, i] += y_series + y_shunt
        yfull[j, j] += y_series + y_shunt
        yfull[i, j] -= y_series
        yfull[j, i] -= y_series
    return yfull


def add_load_admittance(yfull: np.ndarray, pfdata: PFData,
                        include_load: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """把潮流负荷折算为恒阻抗导纳并加到网络。

    使用方法：
        传入 ``Yfull`` 与 ``PFData``，返回 (新导纳矩阵, 每个负荷的 Yload)；
        ``include_load=False`` 时不加负荷、只返回零 Yload。
    步骤：
        每个负荷用 P/V^2 和 -Q/V^2 折算成恒阻抗导纳，累加到对应母线对角。
    """

    result = np.array(yfull, dtype=complex, copy=True)
    yload = np.zeros((pfdata.nload, 3), dtype=float)
    if not include_load:
        return result, yload
    for k, (bus_no, pq, vm) in enumerate(zip(pfdata.load_no, pfdata.load_pq,
                                               pfdata.load_voltage)):
        p_pu, q_pu = pq / pfdata.sbase
        v = vm[0]
        yload[k] = [bus_no, p_pu / (v * v), -q_pu / (v * v)]
        result[int(bus_no) - 1, int(bus_no) - 1] += complex(yload[k, 1], yload[k, 2])
    return result, yload


def kron_reduce(yfull: np.ndarray, gen_buses: Iterable[int],
                active_buses: Optional[Iterable[int]] = None,
                removed_bus: Optional[int] = None) -> Tuple[np.ndarray, ...]:
    """按发电机顺序分块并做 Kron 约简。

    使用方法：
        传入 ``Yfull``、发电机母线号，可选活动母线和被删母线；返回
        (yred, ynn, ynr, yrn, yrr)。
    步骤：
        把发电机节点排到前面、其余排后面，分块后按
        ``Yred = Ynn - Ynr inv(Yrr) Yrn`` 消去负荷节点。
    """

    n = yfull.shape[0]
    if active_buses is None:
        active = np.arange(1, n + 1, dtype=int)
    else:
        active = np.asarray(list(active_buses), dtype=int)
    gen_buses = np.asarray(list(gen_buses), dtype=int)
    gen_pos = [int(np.flatnonzero(active == bus)[0]) for bus in gen_buses]
    load_pos = [i for i in range(n) if i not in gen_pos]
    order = np.asarray(gen_pos + load_pos, dtype=int)
    ordered = yfull[np.ix_(order, order)]
    ngen = len(gen_pos)
    ynn = ordered[:ngen, :ngen]
    ynr = ordered[:ngen, ngen:]
    yrn = ordered[ngen:, :ngen]
    yrr = ordered[ngen:, ngen:]
    if yrr.size:
        yred = ynn - ynr @ np.linalg.solve(yrr, yrn)
    else:
        yred = ynn.copy()
    return yred, ynn, ynr, yrn, yrr


def reorder_structure_preserved(yfull: np.ndarray, gen_buses: Iterable[int],
                                active_buses: Optional[Iterable[int]] = None) -> Tuple[np.ndarray, np.ndarray]:
    """把发电机节点排到前面，供结构保持模型使用。

    使用方法：
        传入 ``Yfull``、发电机母线号和可选活动母线，返回 (重排后的 Yfull,
        Transform)；Transform 记录重排后每个位置对应的原母线编号。
    """

    n = yfull.shape[0]
    active = np.arange(1, n + 1, dtype=int) if active_buses is None else np.asarray(list(active_buses), dtype=int)
    gen_buses = np.asarray(list(gen_buses), dtype=int)
    gen_set = set(gen_buses.tolist())
    non_gen = np.asarray([b for b in active if b not in gen_set], dtype=int)
    transform = np.r_[gen_buses, non_gen]
    index = np.asarray([int(np.flatnonzero(active == b)[0]) for b in transform], dtype=int)
    return yfull[np.ix_(index, index)], transform


def remove_fault_bus(yfull: np.ndarray, active_buses: Iterable[int], fault_bus: int) -> Tuple[np.ndarray, np.ndarray]:
    """删除故障母线所在行列。

    使用方法：
        传入 ``Yfull``、当前活动母线序列和故障母线号，返回 (删后矩阵,
        删后母线序列)。
    """

    active = np.asarray(list(active_buses), dtype=int)
    keep = active != int(fault_bus)
    return yfull[np.ix_(keep, keep)], active[keep]


def remove_fault_line(rxb: np.ndarray, fault_line: Iterable[int]) -> np.ndarray:
    """删除故障切除后应断开的线路（方向不敏感）。

    使用方法：
        传入 RXB 线路矩阵和故障线路 [i, j]，返回删去该线路后的 RXB；正反两个
        方向都能匹配，若找不到该线路则报错。
    """

    line = np.asarray(list(fault_line), dtype=int)
    rxb = np.asarray(rxb, dtype=float)
    mask = ~(
        ((rxb[:, 0].astype(int) == line[0]) & (rxb[:, 1].astype(int) == line[1]))
        | ((rxb[:, 0].astype(int) == line[1]) & (rxb[:, 1].astype(int) == line[0]))
    )
    if mask.all():
        raise ValueError(f"Could not find the line to trip: {line.tolist()}")
    return rxb[mask]
