"""交流潮流与潮流结构转换使用说明。

使用方法：
    调用 ``solve_power_flow(case)`` 求交流潮流，得到 ``PowerFlowResult``；再用
    ``to_pfdata`` 整理成后续网络约简所需的 ``PFData``；如需发电机内电势用
    ``generator_internal_emf``。默认案例不依赖 MATPOWER 或 SciPy。

对应关系：
    本模块承担 MATLAB ``runpf`` 与 ``Fun_ResultBack`` 的项目内等价工作：先解
    PQ/PV/平衡母线潮流，再把电压、发电机、负荷和 RXB 线路数据整理成字段。
"""

from __future__ import annotations

from typing import Tuple

import numpy as np

from .numerics import newton_solve
from .types import CaseData, PFData, PowerFlowResult


def ybus_from_case(case: CaseData) -> np.ndarray:
    """由案例支路参数构造交流潮流用的 Y-bus。

    使用方法：
        传入 ``CaseData``，返回复数导纳矩阵 Y-bus，供潮流残差调用。
    步骤：
        先加母线并联导纳（Gs/Bs 标幺化），再逐条支路加串联导纳、半充电电容和
        变压器变比/移相修正。
    """

    nbus = case.bus.shape[0]
    ybus = np.zeros((nbus, nbus), dtype=complex)
    # 母线 Gs/Bs 的单位是 MW/MVAr at V=1 pu，先换算到标幺导纳再加到对角。
    ybus[np.arange(nbus), np.arange(nbus)] += (
        case.bus[:, 4] + 1j * case.bus[:, 5]
    ) / case.base_mva
    for line in case.branch:
        if line[10] == 0:
            continue
        i, j = int(line[0]) - 1, int(line[1]) - 1
        r, x, b = line[2], line[3], line[4]
        y = 1.0 / complex(r, x)
        tap = line[8] if abs(line[8]) > 1e-12 else 1.0
        shift = np.deg2rad(line[9])
        tap_complex = tap * np.exp(1j * shift)
        ybus[i, i] += (y + 1j * b / 2.0) / (tap * tap)
        ybus[j, j] += y + 1j * b / 2.0
        ybus[i, j] -= y / np.conj(tap_complex)
        ybus[j, i] -= y / tap_complex
    return ybus


def _power_injection(voltage: np.ndarray, ybus: np.ndarray) -> np.ndarray:
    """计算复功率注入 S=V*conj(YV)。

    使用方法：
        传入复电压向量和 Y-bus，返回各母线复功率注入（标幺）。
    """

    return voltage * np.conj(ybus @ voltage)


def solve_power_flow(case: CaseData, tol: float = 1e-10,
                     max_iter: int = 80) -> PowerFlowResult:
    """用有限差分牛顿法求交流潮流。

    使用方法：
        传入 ``CaseData``，返回 ``PowerFlowResult``；``tol`` 控制收敛精度，
        ``max_iter`` 控制最大迭代步数。要求案例恰好一个平衡母线。
    参数：
        case：待求解案例。
        tol：残差收敛阈值。
        max_iter：最大牛顿迭代步数。
    返回：
        含电压、注入、发电机出力、支路潮流和收敛信息的 ``PowerFlowResult``。
    步骤：
        状态量按非平衡母线相角、PQ 母线电压幅值排列；PV 母线只约束有功和电压
        幅值，平衡母线 P/Q 由解反算；收敛后回填发电机出力与支路潮流。
    """

    nbus = case.bus.shape[0]
    bus_type = case.bus[:, 1].astype(int)
    slack = np.flatnonzero(bus_type == 3)
    if slack.size != 1:
        raise ValueError("This teaching power-flow solver requires exactly one slack bus")
    slack = int(slack[0])
    pq = np.flatnonzero(bus_type == 1)
    non_slack = np.array([i for i in range(nbus) if i != slack], dtype=int)
    ybus = ybus_from_case(case)
    p_spec = (np.zeros(nbus) - case.bus[:, 2]) / case.base_mva
    q_spec = (np.zeros(nbus) - case.bus[:, 3]) / case.base_mva
    for row in case.gen:
        bus = int(row[0]) - 1
        p_spec[bus] += row[1] / case.base_mva
        q_spec[bus] += row[2] / case.base_mva
    v0 = case.bus[:, 7] * np.exp(1j * np.deg2rad(case.bus[:, 8]))
    # 发电机电压设定值优先于 bus 表初值。
    for row in case.gen:
        v0[int(row[0]) - 1] = row[5] * np.exp(1j * np.angle(v0[int(row[0]) - 1]))
    x0 = np.r_[np.angle(v0)[non_slack], np.abs(v0)[pq]]

    def unpack(x: np.ndarray) -> np.ndarray:
        # 把状态向量还原成完整复电压：平衡母线固定，PQ 幅值取自状态量。
        angles = np.zeros(nbus)
        magnitudes = np.ones(nbus)
        angles[slack] = np.angle(v0[slack])
        magnitudes[slack] = np.abs(v0[slack])
        angles[non_slack] = x[: non_slack.size]
        magnitudes[non_slack] = np.abs(v0[non_slack])
        magnitudes[pq] = x[non_slack.size :]
        return magnitudes * np.exp(1j * angles)

    def residual(x: np.ndarray) -> np.ndarray:
        # 非平衡母线有功失配 + PQ 母线无功失配组成残差向量。
        inj = _power_injection(unpack(x), ybus)
        p_mis = p_spec[non_slack] - inj.real[non_slack]
        q_mis = q_spec[pq] - inj.imag[pq]
        return np.r_[p_mis, q_mis]

    solution, success, iterations, residual_norm = newton_solve(
        residual, x0, tol=tol, max_iter=max_iter, jac_step=1e-6,
        name="AC power flow"
    )
    voltage = unpack(solution)
    inj = _power_injection(voltage, ybus)
    gen_out = np.zeros_like(case.gen)
    gen_out[:] = case.gen
    for bus in np.unique(case.gen[:, 0].astype(int)):
        # 同一母线多台机组时，按给定出力比例分摊反算的总有功/无功。
        rows = np.flatnonzero(case.gen[:, 0].astype(int) == bus)
        total_p = inj[bus - 1].real * case.base_mva + case.bus[bus - 1, 2]
        total_q = inj[bus - 1].imag * case.base_mva + case.bus[bus - 1, 3]
        p_given = case.gen[rows, 1].sum()
        for r in rows:
            gen_out[r, 1] = total_p * case.gen[r, 1] / p_given if abs(p_given) > 1e-12 else total_p / rows.size
            gen_out[r, 2] = total_q * case.gen[r, 2] / case.gen[rows, 2].sum() if abs(case.gen[rows, 2].sum()) > 1e-12 else total_q / rows.size
    branch_flow = _branch_flow(case, voltage)
    return PowerFlowResult(case, voltage, inj * case.base_mva, gen_out,
                           branch_flow, bool(success), iterations, residual_norm)


def _branch_flow(case: CaseData, voltage: np.ndarray) -> np.ndarray:
    """计算各支路两端潮流。

    使用方法：
        传入案例与已解电压，返回矩阵，列为 [fbus, tbus, Pf, Qf, Pt, Qt]（MW/MVAr）。
    """

    result = np.zeros((case.branch.shape[0], 6), dtype=float)
    for k, line in enumerate(case.branch):
        i, j = int(line[0]) - 1, int(line[1]) - 1
        r, x, b = line[2], line[3], line[4]
        y = 1.0 / complex(r, x)
        tap = line[8] if abs(line[8]) > 1e-12 else 1.0
        shift = np.deg2rad(line[9])
        tapc = tap * np.exp(1j * shift)
        If = (y + 1j * b / 2.0) * voltage[i] / (tap * tap) - y * voltage[j] / np.conj(tapc)
        It = (y + 1j * b / 2.0) * voltage[j] - y * voltage[i] / tapc
        Sf = voltage[i] * np.conj(If) * case.base_mva
        St = voltage[j] * np.conj(It) * case.base_mva
        result[k] = [line[0], line[1], Sf.real, Sf.imag, St.real, St.imag]
    return result


def to_pfdata(result: PowerFlowResult) -> PFData:
    """把潮流结果整理成项目内部 pfdata 结构。

    使用方法：
        传入 ``PowerFlowResult``，返回 ``PFData``；电压整理成 [幅值, 角度] 两列，
        发电机/负荷角度转为 rad，rxb 取支路前 5 列。
    """

    case = result.case
    voltage = np.c_[np.abs(result.voltage), np.rad2deg(np.angle(result.voltage))]
    gen_no = result.gen[:, 0].astype(int)
    gen_pq = result.gen[:, 1:3].copy()
    gen_voltage = voltage[gen_no - 1].copy()
    gen_voltage[:, 1] = np.deg2rad(gen_voltage[:, 1])
    load_rows = np.flatnonzero((case.bus[:, 2] != 0) | (case.bus[:, 3] != 0))
    load_no = case.bus[load_rows, 0].astype(int)
    load_pq = case.bus[load_rows, 2:4].copy()
    load_voltage = voltage[load_rows].copy()
    load_voltage[:, 1] = np.deg2rad(load_voltage[:, 1])
    return PFData(
        sbase=case.base_mva,
        voltage=voltage,
        bus_pq=case.bus[:, 2:4].copy(),
        gen_no=gen_no,
        gen_pq=gen_pq,
        gen_voltage=gen_voltage,
        load_no=load_no,
        load_pq=load_pq,
        load_voltage=load_voltage,
        rxb=case.branch[:, :5].copy(),
        branch_powerflow=result.branch_flow.copy(),
    )


def generator_internal_emf(pfdata: PFData, xd1: np.ndarray,
                           flag_xd: int = 0) -> np.ndarray:
    """计算发电机内电势（对应 Fun_Cal_GenEMF）。

    使用方法：
        传入 ``PFData`` 和暂态电抗 ``xd1``；``flag_xd=0`` 时直接沿用潮流端电压，
        否则由端电压、出力和 xd1 反算内电势幅值与相角。
    返回：
        形状 (ngen, 2) 的 [幅值, 相角(rad)] 矩阵。
    """

    if flag_xd == 0:
        return pfdata.gen_voltage.copy()
    emf = np.zeros((pfdata.ngen, 2), dtype=float)
    for i in range(pfdata.ngen):
        v, angle_deg = pfdata.gen_voltage[i]
        p, q = pfdata.gen_pq[i] / pfdata.sbase
        imag = q * xd1[i] / v
        real = v + imag
        emf[i, 0] = np.sqrt(real * real + (p * xd1[i] / v) ** 2)
        emf[i, 1] = np.arctan2(p * xd1[i] / v, real) + np.deg2rad(angle_deg)
    return emf
