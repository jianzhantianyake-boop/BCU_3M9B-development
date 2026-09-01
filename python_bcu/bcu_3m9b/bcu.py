"""静态初始化与一键实验流程使用说明。

使用方法：
    取默认动态参数用 ``default_preset``；做一次静态初始化（潮流+三工况+SEP）用
    ``build_static_result``，返回 ``StaticResult``；在此基础上跑非交互式 BCU/能量/
    暂态实验用 ``run_bcu_experiment``，返回含退出点、MGP、临界能量、LEA/REA CCT
    的结果字典。
"""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np

from .cases import case9_v2
from .dynamics import find_exitpoint, integrate_reduced, time_domain_cct
from .energy import energy_cct, find_mgp, potential_energy
from .equilibrium import solve_cuep_from_guess, solve_sep
from .network import (add_load_admittance, kron_reduce, remove_fault_bus,
                      remove_fault_line, reorder_structure_preserved,
                      rxb_to_yfull)
from .powerflow import generator_internal_emf, solve_power_flow, to_pfdata
from .types import BaseValue, CaseData, NetworkState, Preset, StaticResult


def default_preset() -> Preset:
    """返回原 Cal_MM_Static 的 3M9B 默认动态参数。

    使用方法：
        直接调用取回 ``Preset``；如需自定义参数，复制返回值再修改对应字段。
    """

    return Preset(
        m=np.array([0.1254, 0.0340, 0.0160]),
        d=np.array([0.01254, 0.0034, 0.0016]),
        pmpu=np.array([0.8980, 1.3432, 0.9419]),
        xd1=np.array([0.0608, 0.1198, 0.1813]),
        epu=np.array([1.1083, 1.1071, 1.0606]),
        path_energy_cal=0,
        equ_cal=2,
        flag_xd=0,
        fault_line=np.array([9, 6], dtype=int),
        fault_position=0,
    )


def _state(name: str, yfull: np.ndarray, active_buses: np.ndarray,
           gen_buses: np.ndarray, removed_bus: Optional[int] = None) -> NetworkState:
    """由完整导纳构造带 Kron 分块的 NetworkState。

    使用方法：
        内部辅助函数；传入工况名、完整导纳、活动母线和发电机母线，返回已做约简
        的 ``NetworkState``。
    """

    yred, ynn, ynr, yrn, yrr = kron_reduce(
        yfull, gen_buses, active_buses=active_buses, removed_bus=removed_bus
    )
    return NetworkState(name, yfull, yred, ynn, ynr, yrn, yrr,
                        removed_bus=removed_bus,
                        metadata={"active_buses": active_buses.copy(), "gen_buses": gen_buses.copy()})


def build_static_result(case: Optional[CaseData] = None,
                        preset: Optional[Preset] = None,
                        solve_tol: float = 1e-9) -> StaticResult:
    """完成对应 MATLAB ``Cal_MM_Static`` 的静态初始化。

    使用方法：
        直接调用取默认 case9_v2 结果，或传入自定义 ``case``、``preset``；返回
        ``StaticResult``。潮流或任一 SEP 不收敛会抛异常。
    参数：
        case：案例数据，缺省用 ``case9_v2()``。
        preset：动态参数，缺省用 ``default_preset()``。
        solve_tol：潮流和 SEP 的收敛阈值。
    步骤：
        (1) 求潮流并整理成 pfdata；(2) 拼预故障导纳并加负荷；(3) 删故障母线得故障
        工况、切线路得故障后工况；(4) 分别求预故障与故障后 SEP；(5) 生成结构保持
        模型所需的发电机优先重排矩阵。
    """

    case = case or case9_v2()
    preset = preset or default_preset()
    basevalue = BaseValue(sbase=case.base_mva)
    pf_result = solve_power_flow(case, tol=solve_tol)
    if not pf_result.success:
        raise RuntimeError(f"Python AC power flow did not converge; residual={pf_result.residual_norm:g}")
    pfdata = to_pfdata(pf_result)
    preset.gen_no = pfdata.gen_no.copy()
    preset.nbus = pfdata.nbus
    preset.s_load = np.c_[pfdata.load_no, pfdata.load_pq / pfdata.sbase]
    emf = generator_internal_emf(pfdata, preset.xd1, preset.flag_xd)
    rxb = pfdata.rxb.copy()
    yorg = rxb_to_yfull(rxb, pfdata.nbus)
    ypre, _ = add_load_admittance(yorg, pfdata)
    active = np.arange(1, pfdata.nbus + 1, dtype=int)
    gen_buses = pfdata.gen_no.copy()
    prefault = _state("prefault", ypre, active, gen_buses)
    prefault.metadata["yorg"] = yorg.copy()
    fault_bus = int(preset.fault_line[preset.fault_position])
    yfault, active_fault = remove_fault_bus(ypre, active, fault_bus)
    fault = _state("fault", yfault, active_fault, gen_buses, removed_bus=fault_bus)
    fault.metadata["yorg"] = yorg[np.ix_(active_fault - 1, active_fault - 1)].copy()
    rxb_post = remove_fault_line(rxb, preset.fault_line)
    ypost_org = rxb_to_yfull(rxb_post, pfdata.nbus)
    ypost, _ = add_load_admittance(ypost_org, pfdata)
    postfault = _state("postfault", ypost, active, gen_buses)
    postfault.metadata["yorg"] = ypost_org.copy()
    delta, omegapu, perr, success, _ = solve_sep(
        preset, prefault, basevalue, np.zeros(preset.ngen), 0.0, tol=solve_tol
    )
    if not success:
        raise RuntimeError(f"Pre-fault SEP did not converge; residual norm={np.linalg.norm(perr):g}")
    delta2, omega2, perr2, success2, _ = solve_sep(
        preset, postfault, basevalue, delta, (omegapu - 1.0) * basevalue.omega_b,
        tol=solve_tol
    )
    if not success2:
        raise RuntimeError(f"Post-fault SEP did not converge; residual norm={np.linalg.norm(perr2):g}")
    # 结构保持模型所需的发电机优先矩阵也在静态初始化阶段一并生成。
    for state in (prefault, fault, postfault):
        modified, transform = reorder_structure_preserved(
            state.yfull, gen_buses, state.metadata["active_buses"]
        )
        yorg_modified, _ = reorder_structure_preserved(
            state.metadata["yorg"], gen_buses, state.metadata["active_buses"]
        )
        state.metadata["yfull_mod"] = modified
        state.metadata["yorg_mod"] = yorg_modified
        state.metadata["transform"] = transform
    return StaticResult(preset, basevalue, pfdata, prefault, fault, postfault,
                        emf, case)


def run_bcu_experiment(static: Optional[StaticResult] = None,
                       fault_time: float = 0.2,
                       tunit: float = 1e-3,
                       postfault_time: float = 2.0,
                       cct_samples: int = 11) -> Dict[str, object]:
    """跑一个不依赖交互输入的 BCU/能量/暂态实验。

    使用方法：
        直接调用会先做静态初始化，或传入已有 ``static``；用 ``fault_time``、
        ``tunit``、``postfault_time``、``cct_samples`` 控制实验设置，返回含轨迹、
        退出点、MGP、CUEP、临界能量和 LEA/REA CCT 的结果字典。
    步骤：
        (1) 以预故障 SEP 为初值积分故障轨迹；(2) 找退出点并追踪 MGP；(3) 以 MGP
        为初值求 CUEP，回落到 SEP 时保留 MGP 作为近似临界点；(4) 算临界能量并用
        能量法得 LEA CCT；(5) 用时域网格得 REA CCT。
    """

    static = static or build_static_result()
    base = static.basevalue
    pre = static.prefault
    fault = static.fault
    post = static.postfault
    delta0 = pre.sep_delta
    omega0 = np.full(static.preset.ngen, pre.sep_omegapu * base.omega_b)
    fault_traj = integrate_reduced(fault_time, tunit, fault, static.preset, base, delta0, omega0)
    exit_index = find_exitpoint(fault_traj, post, static.preset)
    mgp = find_mgp(fault_traj.thetac[exit_index], post, static.preset)  # 保留供诊断/画图/回退
    # CUEP 用 closest-UEP(修正原 find_mgp 求错致 LEA 偏低; 任意 ngen); 失败回退旧 MGP 近似.
    from .cuep import controlling_uep, coi_mismatch
    cres = controlling_uep(static)
    if cres.found:
        critical_delta = cres.cuep
        cuep_source = "closest-UEP (type-1, V>V(SEP))"
        critical_energy = cres.v_cuep
    else:
        critical_delta = mgp["theta_mgp"]
        cuep_source = "MGP fallback (CUEP search failed)"
        critical_energy = float(np.sum(potential_energy(
            static.preset, post, post.sep_delta, critical_delta
        )))
    cuep_perr = coi_mismatch(critical_delta, post.yred, static.preset)
    # LEA 需足够长且较细的故障轨迹以定位能量越界时刻(默认 fault_time 可能短于真实 CCT).
    lea_traj = integrate_reduced(max(fault_time, 0.6), min(tunit, 1e-4),
                                 fault, static.preset, base, delta0, omega0)
    lea = energy_cct(critical_energy, lea_traj, static.preset, post)
    cct_real, real_detail = time_domain_cct(
        fault_time, postfault_time, tunit, fault, post, static.preset,
        base, delta0, omega0, samples=cct_samples
    )
    return {
        "static": static,
        "fault_trajectory": fault_traj,
        "exit_index": exit_index,
        "exit_time": fault_traj.time[exit_index],
        "mgp": mgp,
        "cuep_delta": critical_delta,
        "cuep_perr": cuep_perr,
        "cuep_source": cuep_source,
        "cuep_result": cres,
        "critical_energy": critical_energy,
        "lea": lea,
        "rea_cct": cct_real,
        "rea_detail": real_detail,
    }
