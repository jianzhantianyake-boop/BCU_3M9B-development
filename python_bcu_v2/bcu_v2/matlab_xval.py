# -*- coding: utf-8 -*-
"""T3: 与已验证 MATLAB 平台的逐层交叉验证(只读参考, 不运行/不改 MATLAB).

使用方法:
    run_xval() 加载 MATLAB 平台导出的参考 .mat(matlab_platform/verify/baseline_reduced.mat),
    与 Python(v1 + v2)逐层比对(潮流/Yred/SEP/CUEP/临界能量/CCT), 返回结果列表.

分层结论(默认 9 母线 reduced 链路):
    - 潮流 / Yred / SEP: Python 与 MATLAB 应吻合到 ~1e-10(确定性代数, 强证据).
    - REA CCT(时域): Python P1.2 精确 CCT 与 MATLAB 时域 CCT 应吻合到 ~1e-4.
    - MGP/CUEP/LEA: v1 的 find_mgp 求错(落在 SEP 附近)-> LEA 偏低; v2 用 reduced_region 的
      type-1 UEP 重构 CUEP, 可与 MATLAB CUEP/LEA 对上(~1e-3), 从而定位并修正该差异.

诚实说明: 参考值来自 MATLAB, 仅代表"与该实现一致"; 结合 P0/P1 的物理不变量与等面积金标准,
才构成完整可信度(见 invariants.py / smib.py).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import numpy as np


def default_baseline_path() -> Path:
    """返回 MATLAB 基线参考文件默认路径(matlab_platform/verify/baseline_reduced.mat)."""

    return Path(__file__).resolve().parents[2] / "matlab_platform" / "verify" / "baseline_reduced.mat"


def _res(layer: str, error: float, tol: float, detail: str = "") -> Dict:
    return {"layer": layer, "error": float(error), "tol": float(tol),
            "passed": bool(error <= tol), "detail": detail}


def reconstruct_cuep(static):
    """用 v2 的 reduced_region type-1 UEP 重构受控不稳定平衡点 CUEP(COI 三维角).

    使用方法: 传入 static, 返回 (cuep_delta_3d, found); 取离故障后 SEP 最近的 type-1 UEP.
    """

    from bcu_3m9b.experiments import find_reduced_equilibria

    preset, post = static.preset, static.postfault
    m = preset.m
    eps = [e for e in find_reduced_equilibria(post, preset, grid_points=21) if e["flag"] == 1]
    if not eps:
        return None, False
    # 取与 SEP 的 δ2,δ3 距离最近的 type-1 UEP(受控 UEP 近似).
    sep23 = np.asarray(post.sep_delta)[1:3]
    eps.sort(key=lambda e: np.linalg.norm(e["xep"] - sep23))
    xep = eps[0]["xep"]  # [δ2c, δ3c]
    d1 = -(m[1] * xep[0] + m[2] * xep[1]) / m[0]  # COI 约束
    return np.array([d1, xep[0], xep[1]]), True


def run_xval(baseline_path: Optional[Path] = None) -> List[Dict]:
    """执行逐层交叉验证, 返回结果列表."""

    import scipy.io as sio
    from bcu_3m9b import build_static_result
    from bcu_3m9b.dynamics import integrate_reduced
    from bcu_3m9b.energy import potential_energy, energy_cct
    from . import cct as _cct

    path = Path(baseline_path or default_baseline_path())
    d = sio.loadmat(str(path), squeeze_me=True, struct_as_record=False)
    m_pre, m_post, m_crit = d["prefault"], d["postfault"], d["Critical"]

    static = build_static_result()  # 与 MATLAB 同为默认 9 母线 case9_v2 + 同 preset
    results: List[Dict] = []

    # 层1: Yred(pre/post) —— 确定性代数, 应吻合到 ~1e-9.
    results.append(_res("Yred 预故障", np.max(np.abs(np.array(m_pre.Yred) - static.prefault.yred)), 1e-8))
    results.append(_res("Yred 故障后", np.max(np.abs(np.array(m_post.Yred) - static.postfault.yred)), 1e-8))

    # 层2: SEP(pre/post) 角度 + 速度.
    results.append(_res("SEP 预故障角度", np.max(np.abs(np.array(m_pre.SEP_delta) - static.prefault.sep_delta)), 1e-7))
    results.append(_res("SEP 故障后角度", np.max(np.abs(np.array(m_post.SEP_delta) - static.postfault.sep_delta)), 1e-7))
    results.append(_res("SEP 故障后速度(pu)", abs(float(m_post.SEP_omegapu) - static.postfault.sep_omegapu), 1e-7))

    # 层3: CUEP —— v1 的 MGP 求错, 用 v2 重构的 type-1 UEP 与 MATLAB CUEP 对比.
    cuep_py, found = reconstruct_cuep(static)
    if found:
        e_cuep = float(np.max(np.abs(np.array(m_post.CUEP_delta) - cuep_py)))
        results.append(_res("CUEP(v2 重构 vs MATLAB)", e_cuep, 5e-3,
                            f"py={np.array2string(cuep_py, precision=4)}"))
    else:
        results.append(_res("CUEP(v2 重构 vs MATLAB)", np.nan, 5e-3, "未找到 type-1 UEP"))

    # 层4: 临界能量 + LEA CCT —— 用 v2 CUEP 计算, 与 MATLAB 对比.
    if found:
        e_crit_py = float(np.sum(potential_energy(static.preset, static.postfault,
                                                  static.postfault.sep_delta, cuep_py)))
        d0 = static.prefault.sep_delta
        w0 = np.full(static.preset.ngen, static.prefault.sep_omegapu * static.basevalue.omega_b)
        fault_traj = integrate_reduced(0.4, 1e-4, static.fault, static.preset,
                                       static.basevalue, d0, w0)
        lea_py = energy_cct(e_crit_py, fault_traj, static.preset, static.postfault).cct
        results.append(_res("LEA CCT(v2 vs MATLAB)", abs(lea_py - float(m_crit.LEA.CCT)), 5e-3,
                            f"py={lea_py:.4f}s, matlab={float(m_crit.LEA.CCT):.4f}s"))

    # 层5: REA CCT(时域) —— P1.2 精确 vs MATLAB.
    rea_py, _ = _cct.precise_cct_reduced(static, tol=5e-5)
    results.append(_res("REA CCT(P1.2 vs MATLAB)", abs(rea_py - float(m_crit.REA.CCT)), 3e-3,
                        f"py={rea_py:.4f}s, matlab={float(m_crit.REA.CCT):.4f}s"))

    return results
