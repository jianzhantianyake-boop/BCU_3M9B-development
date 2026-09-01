# -*- coding: utf-8 -*-
"""通用能量法 CUEP + LEA CCT(任意 ngen 的 closest-UEP 法).

背景:
    v2 之前只有 3 机专用的 CUEP 重构(matlab_xval.reconstruct_cuep, 靠 experiments 的二维网格
    枚举 type-1 UEP, 取离 SEP 最近的). 该网格法在 39 母线(10 机, 9 维)会组合爆炸, 于是 39 母线
    出不了能量法 LEA. 本模块把 CUEP 的确定推广到**任意机数**.

方法(closest controlling UEP, 与 MATLAB 平台一致):
    1. 从**结构化初值**做 Newton 求梯度系统零点 F(θ)=0(F 为 COI 功率失配 = 故障后势能负梯度):
       - MOD(Mode of Disturbance)翻转: 对每个不超过 max_group 台机的子集 S, 把 S 内各机相角
         镜像/平移约 π(θ_j -> π-θ_j 或 θ_j+π), 其余保持 SEP -> 覆盖"一群机相对其余机失步"的
         各阶不稳定模式;
       - 故障轨迹采样: 取 fault-on 轨迹上若干点作补充初值, 捕捉与该故障相关的 UEP.
    2. 每个收敛解校验为 type-1(去平凡零方向后恰一个正特征值)且 V(θ)>V(SEP), 去重后得 type-1
       UEP 集合;
    3. 取**离故障后 SEP 最近**的 type-1 UEP 作为 CUEP(closest-UEP 准则). 其临界势能 V(CUEP)-V(SEP)
       即能量法临界能量.

为何不用 BCU exit-point 法:
    实测本平台默认故障([9,6])下, fault-on 轨迹的 PEBS 退出点判据退化(功率失配·相对速度在故障时
    长内仅在 0.34s 处过零, 远晚于 CCT≈0.24s), 且 fault-on 轨迹不指向该 closest UEP -> 从退出点
    正向积分梯度系统被 CUEP 的不稳定方向推飞. 而 MATLAB 平台本身用的就是 closest-UEP(网格枚举取
    最近), 故通用化沿 closest-UEP 路线, 既与 T3 一致又可扩展. exit-point 法留作后续严格 BCU 的
    方向(见交接文档第 8 节).

正确性锚点:
    - 3 机: 本模块 CUEP 与 MATLAB CUEP 对到 ~1e-11, LEA 0.2274 vs 0.2275(见 test_cuep.py);
    - 多机: 以物理不变量自检 —— type-1(恰一个正特征值) + V(CUEP)>V(SEP) + LEA<=REA(时域真值).
    诚实说明: closest-UEP 是 controlling UEP 的常用工程近似, 不保证在强非线性/多摆工况总等于严格
    controlling UEP; 定量结论请结合时域 REA(cct.precise_cct_reduced)交叉核对.

单位: 角度 rad, 时间 s, 功率/导纳 pu. 全程 COI 坐标(Σ m_i θ_i = 0).
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import List, Optional

import numpy as np


# ------------------------------- 梯度系统 -------------------------------

def coi_mismatch(theta: np.ndarray, yred: np.ndarray, preset) -> np.ndarray:
    """通用 COI 功率失配 F(θ)(梯度系统右端, 零点即平衡点, 任意 ngen).

    使用方法:
        传入功角 θ(rad, 长度 n)、约简导纳 yred、参数 preset, 返回 n 维
        F_i = Pm_i - Pe_i(θ) - (m_i/M)·Pcoi, Pcoi = Σ(Pm - Pe), M = Σ m.
    说明:
        F 是故障后势能在 COI 度量下的负梯度; F(θ)=0 即网络约简模型的平衡点(SEP / type-k UEP).
    """

    from bcu_3m9b.equilibrium import electrical_power

    theta = np.asarray(theta, dtype=float).reshape(-1)
    m = np.asarray(preset.m, dtype=float)
    pe = electrical_power(theta, yred, preset.epu)
    pcoi = np.sum(preset.pmpu - pe)
    return preset.pmpu - pe - (m / np.sum(m)) * pcoi


def _project_coi(theta: np.ndarray, m: np.ndarray) -> np.ndarray:
    """把角度投影到 COI 流形 Σ m_i θ_i = 0(减去质量加权平均)."""

    return theta - np.dot(m, theta) / np.sum(m)


def _reduced_jacobian_eig(theta: np.ndarray, yred: np.ndarray, preset) -> np.ndarray:
    """在 COI 子空间上求梯度系统 Jacobian 的特征值(实部, 去平凡零方向).

    使用方法:
        返回 n-1 个特征值实部, 用于判定平衡点类型. 平凡零方向来自 electrical_power 只依赖角度差
        (整体平移 θ+c·1 不改 F), 对应一个恒零特征值, 按 |实部| 最小剔除.
    """

    from bcu_3m9b.numerics import numerical_jacobian

    jac = numerical_jacobian(lambda z: coi_mismatch(z, yred, preset), theta)
    lam = np.linalg.eigvals(jac)
    order = np.argsort(np.abs(np.real(lam)))
    return np.real(lam[order[1:]])


# --------------------------- type-1 UEP 搜索 ---------------------------

def _newton_equilibrium(static, guess: np.ndarray, tol: float = 1e-10):
    """以 guess 为初值无副作用地求平衡点; 返回 (theta_coi, resid) 或 None.

    用 scipy.optimize.root(hybr, 见 solvers.solve_sep_scipy)加速: 比 v1 数值雅可比牛顿快约一个
    数量级, 使多机(上百初值)搜索可行, 且无副作用(不污染 post.SEP).
    """

    from .solvers import solve_sep_scipy

    preset, base, post = static.preset, static.basevalue, static.postfault
    sol, _, _, ok = solve_sep_scipy(preset, post, base, delta0=np.asarray(guess, float),
                                    omega0=0.0, tol=tol)
    sol = _project_coi(np.asarray(sol, dtype=float), preset.m)
    resid = float(np.linalg.norm(coi_mismatch(sol, post.yred, preset)))
    if not (ok and resid < 1e-7):
        return None
    return sol, resid


def find_type1_ueps(static, max_group: int = 2, fault_samples: int = 8,
                    tfault: float = 0.6, tunit: float = 1e-4,
                    dedup_tol: float = 1e-3) -> List[dict]:
    """搜索故障后系统的 type-1 UEP 集合(结构化初值 Newton, 任意 ngen).

    使用方法:
        传入 v1 的 StaticResult, 返回 type-1 UEP 列表(按离 SEP 距离升序), 每项字典含:
        theta(COI 角), v(V(θ)-V(SEP)), dist(到 SEP 的欧氏距离), eig(约简特征值实部).
    参数:
        max_group: MOD 初值中联合翻转的最大机数(controlling UEP 通常为低阶模式, 3 足够).
        fault_samples: 额外从 fault-on 轨迹采样的初值个数(捕捉故障相关 UEP).
        tfault/tunit: 故障轨迹时长与步长(仅用于生成采样初值).
    步骤:
        (1) 生成初值: MOD 子集翻转(mirror: θ_j->π-θ_j; plus: θ_j+π) + 故障轨迹采样;
        (2) 每个初值 Newton 求 F=0; (3) 校验 type-1 且 V>V(SEP), 去重.
    """

    from bcu_3m9b.dynamics import integrate_reduced
    from bcu_3m9b.energy import potential_energy

    preset, base, post = static.preset, static.basevalue, static.postfault
    n = preset.ngen
    sep = np.asarray(post.sep_delta, dtype=float)

    guesses: List[np.ndarray] = []
    # (1a) MOD 翻转初值.
    idx = list(range(n))
    for r in range(1, min(max_group, n) + 1):
        for S in itertools.combinations(idx, r):
            for mode in ("mirror", "plus"):
                g = sep.copy()
                for j in S:
                    g[j] = (np.pi - sep[j]) if mode == "mirror" else (sep[j] + np.pi)
                guesses.append(_project_coi(g, preset.m))
    # (1b) 故障轨迹采样初值.
    if fault_samples > 0:
        d0 = np.asarray(static.prefault.sep_delta, dtype=float)
        w0 = np.full(n, static.prefault.sep_omegapu * base.omega_b)
        ftraj = integrate_reduced(tfault, tunit, static.fault, preset, base, d0, w0)
        for k in np.linspace(0, ftraj.time.size - 1, fault_samples).astype(int):
            guesses.append(_project_coi(ftraj.theta[k], preset.m))

    found: List[dict] = []
    for g in guesses:
        out = _newton_equilibrium(static, g)
        if out is None:
            continue
        sol, _ = out
        if np.linalg.norm(sol - sep) < 1e-3:
            continue  # 回落到 SEP
        eigr = _reduced_jacobian_eig(sol, post.yred, preset)
        if int(np.sum(eigr > 1e-6)) != 1:
            continue  # 只保留 type-1
        v = float(np.sum(potential_energy(preset, post, sep, sol)))
        if v <= 0.0:
            continue  # 要求 V(UEP) > V(SEP)
        if any(np.max(np.abs(sol - f["theta"])) < dedup_tol for f in found):
            continue
        found.append({"theta": sol, "v": v, "dist": float(np.linalg.norm(sol - sep)),
                      "eig": eigr})
    found.sort(key=lambda f: f["dist"])
    return found


# ------------------------------- CUEP -------------------------------

@dataclass
class CUEPResult:
    """CUEP 求解结果.

    字段:
        cuep: 受控 UEP 的 COI 角(n 维, rad). found: 是否成功.
        v_cuep: 临界势能 V(CUEP)-V(SEP). dist: CUEP 到 SEP 的距离. eig_reduced: 约简特征值实部.
        n_type1: 找到的 type-1 UEP 总数. note: 诊断.
    """

    cuep: Optional[np.ndarray]
    found: bool
    v_cuep: float = float("nan")
    dist: float = float("nan")
    eig_reduced: Optional[np.ndarray] = None
    n_type1: int = 0
    note: str = ""


def controlling_uep(static, max_group: int = 2, fault_samples: int = 8) -> CUEPResult:
    """求任意 ngen 的受控不稳定平衡点 CUEP(closest-UEP 准则).

    使用方法:
        传入 v1 的 StaticResult, 返回 CUEPResult; 取 find_type1_ueps 里离 SEP 最近的 type-1 UEP.
    """

    ueps = find_type1_ueps(static, max_group=max_group, fault_samples=fault_samples)
    if not ueps:
        return CUEPResult(None, False, note="no type-1 UEP with V>V(SEP) found")
    best = ueps[0]
    return CUEPResult(best["theta"], True, v_cuep=best["v"], dist=best["dist"],
                      eig_reduced=best["eig"], n_type1=len(ueps),
                      note="closest type-1 UEP (controlling UEP approximation)")


# ------------------------------- LEA CCT -------------------------------

@dataclass
class LEAResult:
    """能量法 LEA CCT 结果.

    字段: lea(s), found, critical_energy(V(CUEP)-V(SEP)), cuep(CUEPResult), note.
    """

    lea: float
    found: bool
    critical_energy: float = float("nan")
    cuep: Optional[CUEPResult] = None
    note: str = ""


def energy_lea_cct(static, tfault: float = 0.6, tunit: float = 1e-4,
                   max_group: int = 2, fault_samples: int = 8) -> LEAResult:
    """任意 ngen 的能量法 LEA CCT.

    使用方法:
        传入 StaticResult, 返回 LEAResult. 先用 controlling_uep 求 CUEP 得临界能量, 再沿故障轨迹
        用 energy_cct 找总能量首次越过临界能量的时刻 = LEA CCT.
    参数:
        tfault/tunit: 故障轨迹时长与步长(需 tfault > 预计 CCT).
    """

    from bcu_3m9b.dynamics import integrate_reduced
    from bcu_3m9b.energy import energy_cct

    cres = controlling_uep(static, max_group=max_group, fault_samples=fault_samples)
    if not cres.found or cres.cuep is None:
        return LEAResult(float("nan"), False, cuep=cres, note=f"CUEP failed: {cres.note}")

    preset, base = static.preset, static.basevalue
    ecrit = cres.v_cuep
    d0 = np.asarray(static.prefault.sep_delta, dtype=float)
    w0 = np.full(preset.ngen, static.prefault.sep_omegapu * base.omega_b)
    fault_traj = integrate_reduced(tfault, tunit, static.fault, preset, base, d0, w0)
    lea = energy_cct(ecrit, fault_traj, preset, static.postfault)
    return LEAResult(float(lea.cct), bool(lea.flag_cct), critical_energy=ecrit,
                     cuep=cres, note="" if lea.flag_cct else "energy never crossed critical")
