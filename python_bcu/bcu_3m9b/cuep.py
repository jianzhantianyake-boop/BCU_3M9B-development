# -*- coding: utf-8 -*-
"""通用能量法 CUEP + LEA CCT(任意 ngen 的 closest-UEP 法)。

背景：
    修正原 v1 隐患——``energy.find_mgp`` 的最小梯度点追踪步长/步数不足, MGP 落在 SEP 附近,
    致临界能量偏低、能量法 LEA CCT 偏低。本模块用 closest-UEP 法直接确定受控不稳定平衡点 CUEP,
    对任意机数成立(不再假设 ngen=3)。

方法(closest controlling UEP)：
    (1) 从结构化初值做求根 F(θ)=0(F 为 COI 功率失配 = 故障后势能负梯度)：
        - MOD(Mode of Disturbance)翻转: 对不超过 max_group 台机的子集 S, 把 S 内各机相角镜像/
          平移约 π(θ_j -> π-θ_j 或 θ_j+π), 其余保持 SEP, 覆盖"一群机相对其余机失步"的各阶模式;
        - 故障轨迹采样: 取 fault-on 轨迹若干点作补充初值, 捕捉故障相关 UEP.
    (2) 每个解校验为 type-1(去平凡零方向后恰一个正特征值)且 V(θ)>V(SEP), 去重;
    (3) 取离故障后 SEP 最近的 type-1 UEP 作为 CUEP。其临界势能 V(CUEP)-V(SEP) 即临界能量。

诚实说明：
    closest-UEP 是 controlling UEP 的常用工程近似, 不保证在强非线性/多摆工况总等于严格
    controlling UEP; 定量结论请结合时域真值(dynamics.time_domain_cct)交叉核对。9 母线上本法
    与 MATLAB 平台 CUEP 一致到 ~1e-11。

单位: 角度 rad, 时间 s, 功率/导纳 pu。全程 COI 坐标(Σ m_i θ_i = 0)。依赖 scipy。
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import List, Optional

import numpy as np

from .equilibrium import electrical_power, normalize_coi, sep_residual
from .numerics import numerical_jacobian


# ------------------------------- 梯度系统 -------------------------------

def coi_mismatch(theta: np.ndarray, yred: np.ndarray, preset) -> np.ndarray:
    """通用 COI 功率失配 F(θ)=Pm-Pe(θ)-(m/M)·Pcoi(梯度系统右端, 零点即平衡点, 任意 ngen)。"""

    theta = np.asarray(theta, dtype=float).reshape(-1)
    m = np.asarray(preset.m, dtype=float)
    pe = electrical_power(theta, yred, preset.epu)
    pcoi = np.sum(preset.pmpu - pe)
    return preset.pmpu - pe - (m / np.sum(m)) * pcoi


def _project_coi(theta: np.ndarray, m: np.ndarray) -> np.ndarray:
    """把角度投影到 COI 流形 Σ m_i θ_i = 0。"""

    return theta - np.dot(m, theta) / np.sum(m)


def _reduced_jacobian_eig(theta: np.ndarray, yred: np.ndarray, preset) -> np.ndarray:
    """COI 子空间上梯度系统 Jacobian 特征值(实部, 去平凡零方向); 正特征值个数=不稳定维数。"""

    jac = numerical_jacobian(lambda z: coi_mismatch(z, yred, preset), theta)
    lam = np.linalg.eigvals(jac)
    order = np.argsort(np.abs(np.real(lam)))
    return np.real(lam[order[1:]])


def _solve_equilibrium(preset, state, basevalue, guess, tol: float = 1e-10):
    """以 guess 为初值 scipy.root 求平衡点(无副作用); 返回 (theta_coi, resid) 或 None。"""

    from scipy.optimize import root

    n = preset.ngen
    g = normalize_coi(np.asarray(guess, dtype=float), preset.m)
    z0 = np.r_[g[: n - 1] - g[-1], 0.0]
    sol = root(lambda z: sep_residual(z, preset, state, basevalue), z0, method="hybr", tol=tol)
    theta = normalize_coi(np.r_[sol.x[: n - 1], 0.0], preset.m)
    resid = float(np.linalg.norm(coi_mismatch(theta, state.yred, preset)))
    if not (sol.success and resid < 1e-7):
        return None
    return theta, resid


# --------------------------- type-1 UEP 搜索 ---------------------------

def find_type1_ueps(static, max_group: int = 2, fault_samples: int = 8,
                    tfault: float = 0.6, tunit: float = 1e-4,
                    dedup_tol: float = 1e-3) -> List[dict]:
    """搜索故障后系统的 type-1 UEP 集合(结构化初值求根, 任意 ngen)。

    返回按离 SEP 距离升序的列表, 每项字典含 theta(COI 角)、v(V-V(SEP))、dist、eig(约简特征值)。
    """

    from .dynamics import integrate_reduced
    from .energy import potential_energy

    preset, base, post = static.preset, static.basevalue, static.postfault
    n = preset.ngen
    sep = np.asarray(post.sep_delta, dtype=float)

    guesses: List[np.ndarray] = []
    idx = list(range(n))
    for r in range(1, min(max_group, n) + 1):
        for S in itertools.combinations(idx, r):
            for mode in ("mirror", "plus"):
                g = sep.copy()
                for j in S:
                    g[j] = (np.pi - sep[j]) if mode == "mirror" else (sep[j] + np.pi)
                guesses.append(_project_coi(g, preset.m))
    if fault_samples > 0:
        d0 = np.asarray(static.prefault.sep_delta, dtype=float)
        w0 = np.full(n, static.prefault.sep_omegapu * base.omega_b)
        ftraj = integrate_reduced(tfault, tunit, static.fault, preset, base, d0, w0)
        for k in np.linspace(0, ftraj.time.size - 1, fault_samples).astype(int):
            guesses.append(_project_coi(ftraj.theta[k], preset.m))

    found: List[dict] = []
    for g in guesses:
        out = _solve_equilibrium(preset, post, base, g)
        if out is None:
            continue
        sol, _ = out
        if np.linalg.norm(sol - sep) < 1e-3:
            continue
        eigr = _reduced_jacobian_eig(sol, post.yred, preset)
        if int(np.sum(eigr > 1e-6)) != 1:
            continue
        v = float(np.sum(potential_energy(preset, post, sep, sol)))
        if v <= 0.0:
            continue
        if any(np.max(np.abs(sol - f["theta"])) < dedup_tol for f in found):
            continue
        found.append({"theta": sol, "v": v, "dist": float(np.linalg.norm(sol - sep)),
                      "eig": eigr})
    found.sort(key=lambda f: f["dist"])
    return found


# ------------------------------- CUEP -------------------------------

@dataclass
class CUEPResult:
    """CUEP 求解结果: cuep(COI 角), found, v_cuep(V(CUEP)-V(SEP)), dist, eig_reduced, n_type1, note。"""

    cuep: Optional[np.ndarray]
    found: bool
    v_cuep: float = float("nan")
    dist: float = float("nan")
    eig_reduced: Optional[np.ndarray] = None
    n_type1: int = 0
    note: str = ""


def controlling_uep(static, max_group: int = 2, fault_samples: int = 8) -> CUEPResult:
    """求任意 ngen 的受控不稳定平衡点 CUEP(取离 SEP 最近的 type-1 UEP)。"""

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
    """能量法 LEA CCT 结果: lea(s), found, critical_energy(V(CUEP)-V(SEP)), cuep(CUEPResult), note。"""

    lea: float
    found: bool
    critical_energy: float = float("nan")
    cuep: Optional[CUEPResult] = None
    note: str = ""


def energy_lea_cct(static, tfault: float = 0.6, tunit: float = 1e-4,
                   max_group: int = 2, fault_samples: int = 8) -> LEAResult:
    """任意 ngen 的能量法 LEA CCT: 求 CUEP 得临界能量, 再沿故障轨迹找总能量首次越界时刻。"""

    from .dynamics import integrate_reduced
    from .energy import energy_cct

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
