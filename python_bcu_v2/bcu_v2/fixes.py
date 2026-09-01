# -*- coding: utf-8 -*-
"""P0.2: 修正 v1 的 5 个已知隐患(不改冻结的 v1, 在 v2 提供正确行为).

隐患与对策:
    1. solve_sep 有副作用(改写 NetworkState.sep_*)  -> solve_sep_pure: 存/复原, 不留痕.
    2. run_bcu_experiment 污染输入 static 的 postfault.SEP -> run_experiment_clean: 深拷贝再算.
    3. find_exitpoint 初始伪过零(退出点 index=0)       -> find_exitpoint_fixed: 排除初始抖动.
    4. trajectory_stable 轻阻尼判据失效                 -> is_stable_bounded: 有界性判据.
    5. SPM 代数解偶发不收敛                             -> solve_algebraic_robust: scipy.root.
"""

from __future__ import annotations

import copy
from typing import Tuple

import numpy as np


def solve_sep_pure(preset, state, basevalue, **kwargs):
    """无副作用地求 SEP: 调用 v1.solve_sep 后复原 state 的 sep_* 字段.

    使用方法: 与 v1.solve_sep 同参, 返回 (delta, omegapu, perr, success, iters), 但**不**改写
    传入的 state(适合在扫描/共用同一 state 时避免污染).
    """

    from bcu_3m9b.equilibrium import solve_sep

    saved = (state.sep_delta, state.sep_omegapu, state.sep_perr)
    try:
        return solve_sep(preset, state, basevalue, **kwargs)
    finally:
        state.sep_delta, state.sep_omegapu, state.sep_perr = saved


def run_experiment_clean(static=None, **kwargs):
    """跑 BCU 实验但不污染传入的 static(在深拷贝上运行).

    使用方法: 缺省会新建 static; 传入的 static 保持不变(其 postfault.SEP 不会被 CUEP 改写).
    返回 v1.run_bcu_experiment 的结果字典.
    """

    from bcu_3m9b.bcu import build_static_result, run_bcu_experiment

    static = static or build_static_result()
    return run_bcu_experiment(copy.deepcopy(static), **kwargs)


def find_exitpoint_fixed(traj, postfault, preset, eps: float = 1e-9) -> int:
    """修正版退出点: 排除初始 ωc≈0 造成的伪过零.

    使用方法: 与 v1.find_exitpoint 同参, 返回退出点索引; 要求过零点处功率失配与相对速度的
    点积由"明显为负"变"明显为正", 且跳过第 0 个点.
    """

    from bcu_3m9b.equilibrium import electrical_power

    values = np.zeros(traj.time.size)
    for k in range(traj.time.size):
        pe = electrical_power(traj.theta[k], postfault.yred, preset.epu)
        values[k] = np.dot(preset.pmpu - pe, traj.omegac[k])
    mask = (values[:-1] < -eps) & (values[1:] > eps)
    if mask.size:
        mask[0] = False  # 排除初始伪过零
    cross = np.flatnonzero(mask)
    return int(cross[0]) if cross.size else traj.time.size - 1


def is_stable_bounded(traj, sep_delta, limit: float = 1.5 * np.pi) -> bool:
    """有界性稳定判据(适合轻阻尼): 全程有限且离 SEP 不超过 limit 即稳定.

    使用方法: 传入轨迹与故障后 SEP 角, 返回 True/False. 替代 v1 对轻阻尼失效的"回到 SEP".
    """

    dc = np.asarray(traj.thetac)
    return bool(np.all(np.isfinite(dc)) and np.max(np.abs(dc - sep_delta)) < limit)


def solve_algebraic_robust(delta_gen, state, preset, guess=None,
                           tol: float = 1e-10) -> Tuple[np.ndarray, bool, float]:
    """用 scipy.optimize.root 稳健求解 SPM 负荷代数方程(替 v1 自写牛顿).

    使用方法: 与 v1.solve_algebraic 同接口, 返回 (解, 是否收敛, 残差范数); 对 v1 偶发首步
    不收敛的工况更稳健.
    """

    from scipy.optimize import root
    from bcu_3m9b.spm import algebraic_residual, _load_power

    yorg = np.asarray(state.metadata.get("yorg_mod", state.yfull), dtype=complex)
    transform = np.asarray(state.metadata.get("transform"), dtype=int)
    ngen = preset.ngen
    load_pq = _load_power(preset, transform[ngen:])
    nload = yorg.shape[0] - ngen
    if guess is None:
        guess = np.r_[np.zeros(nload), np.ones(nload)]

    def resid(z):
        return algebraic_residual(z, delta_gen, yorg, load_pq, ngen)

    sol = root(resid, guess, method="hybr", tol=tol)
    res_norm = float(np.linalg.norm(resid(sol.x)))
    return sol.x, bool(sol.success and res_norm < 1e-6), res_norm
