# -*- coding: utf-8 -*-
"""P0.1: 不依赖 MATLAB 的正确性校验(T1 物理不变量 + T2 独立参照).

使用方法:
    调用 run_all() 跑全部检查, 返回结果列表, 每项 dict 含 name/passed/error/tol/detail.
    run_validation.py 会把它渲染成表格.

各检查(可信度从高到低):
    T1: 潮流残差, SEP 平衡与稳定, 能量单调, CCT 夹逼(LEA<=真值), 极限工况.
    T2: SMIB 等面积闭式 CCT vs 数值 CCT(金标准), 轨迹 vs scipy 积分器交叉.
"""

from __future__ import annotations

from typing import Dict, List

import numpy as np

from . import smib as _smib
from . import cct as _cct
from . import fixes as _fixes


def _result(name: str, passed: bool, error: float, tol: float, detail: str = "") -> Dict:
    return {"name": name, "passed": bool(passed), "error": float(error),
            "tol": float(tol), "detail": detail}


# ------------------------------- T1 -------------------------------

def check_powerflow_residual(tol: float = 1e-6) -> Dict:
    """T1: 交流潮流功率失配应趋于 0."""

    from bcu_3m9b import case9_v2, solve_power_flow
    pf = solve_power_flow(case9_v2(), tol=1e-9)
    return _result("潮流功率失配 -> 0", pf.residual_norm < tol, pf.residual_norm, tol,
                   f"iters={pf.iterations}")


def check_sep_equilibrium(static, tol: float = 1e-6) -> Dict:
    """T1: 预/故障后 SEP 的功率残差应趋于 0(是真平衡点)."""

    e = max(np.linalg.norm(static.prefault.sep_perr),
            np.linalg.norm(static.postfault.sep_perr))
    return _result("SEP 是真平衡点(残差->0)", e < tol, e, tol)


def check_sep_stable(static, tol: float = 0.0) -> Dict:
    """T1: 故障后 SEP 处梯度系统雅可比特征值应全部 < 0(稳定结点)."""

    from bcu_3m9b.experiments import reduced_gradient
    from bcu_3m9b.numerics import numerical_jacobian
    post, preset = static.postfault, static.preset
    deltac = np.asarray(post.sep_delta)[1:3]
    jac = numerical_jacobian(lambda z: reduced_gradient(z, post, preset), deltac)
    max_re = float(np.max(np.real(np.linalg.eigvals(jac))))
    return _result("SEP 是稳定结点(Re(λ)<0)", max_re < 0.0, max_re, 0.0,
                   f"max Re(λ)={max_re:.3g}")


def check_energy_monotonic(static, tol: float = 1e-3) -> Dict:
    """T1: 有阻尼下故障后总能量应单调不增(允许微小数值抖动)."""

    from bcu_3m9b.dynamics import integrate_reduced
    from bcu_3m9b.energy import trajectory_energy
    preset, base = static.preset, static.basevalue
    d0 = static.prefault.sep_delta  # 相对故障后 SEP 有初始势能
    w0 = np.full(preset.ngen, static.postfault.sep_omegapu * base.omega_b)
    traj = integrate_reduced(1.0, 1e-3, static.postfault, preset, base, d0, w0)
    total = trajectory_energy(traj, preset, static.postfault)["total"]
    scale = max(1e-9, abs(total[0]))
    max_rise = float(np.max(np.diff(total))) / scale  # 相对最大上升
    net = (total[-1] - total[0]) / scale
    return _result("总能量单调不增(阻尼耗散)", max_rise < tol and net <= tol, max_rise, tol,
                   f"净变化={net:.2e}(应<=0)")


def check_cct_sandwich(static, tol: float = 5e-3) -> Dict:
    """T1: 能量法 LEA CCT 应是保守下界, 即 LEA <= 精确时域 CCT."""

    lea = _fixes.run_experiment_clean(static)["lea"].cct
    rea, found = _cct.precise_cct_reduced(static)
    return _result("CCT 夹逼: LEA <= 精确REA", lea <= rea + tol, lea - rea, tol,
                   f"LEA={lea:.4g}s, 精确REA={rea:.4g}s, found={found}")


def check_limit_cases(static) -> Dict:
    """T1: 零故障必稳, 极长故障必失稳."""

    from bcu_3m9b.dynamics import integrate_reduced
    preset, base = static.preset, static.basevalue
    d0 = static.prefault.sep_delta
    w0 = np.full(preset.ngen, static.prefault.sep_omegapu * base.omega_b)

    def stable_after(tc):
        f = integrate_reduced(max(tc, 1e-6), 1e-3, static.fault, preset, base, d0, w0)
        p = integrate_reduced(2.0, 1e-3, static.postfault, preset, base, f.theta[-1], f.omega[-1])
        return _fixes.is_stable_bounded(p, static.postfault.sep_delta)

    zero_ok = stable_after(0.0)
    long_ok = not stable_after(1.0)
    return _result("极限: 零故障稳 & 长故障失稳", zero_ok and long_ok,
                   0.0 if (zero_ok and long_ok) else 1.0, 0.0,
                   f"零故障稳={zero_ok}, 长故障失稳={long_ok}")


# ------------------------------- T2 -------------------------------

def check_smib_analytic(rel_tol: float = 1e-2) -> Dict:
    """T2(金标准): SMIB 数值精确 CCT vs 等面积闭式 CCT(母线金属性短路).

    同时校验数值切除角 ≈ 闭式临界切除角. 相对误差应 < rel_tol.
    """

    from scipy.integrate import solve_ivp
    s = _smib.default_smib(bolted=True)
    tcc_analytic = _smib.critical_time_analytic(s)
    dcc_analytic = _smib.critical_angle(s)
    tcc_numeric, found = _cct.precise_cct_smib(s, tol=1e-6)
    # 数值 CCT 时刻对应的切除角
    d0 = _smib.delta0(s)
    sol = solve_ivp(lambda t, z: _smib.swing_rhs(z, s, "fault"), [0, tcc_numeric],
                    [d0, 0.0], rtol=1e-10, atol=1e-12)
    dcc_numeric = float(sol.y[0, -1])
    err_t = abs(tcc_numeric - tcc_analytic) / abs(tcc_analytic)
    err_d = abs(dcc_numeric - dcc_analytic) / abs(dcc_analytic)
    err = max(err_t, err_d)
    return _result("SMIB: 数值CCT vs 等面积闭式", found and err < rel_tol, err, rel_tol,
                   f"tcc: 数值={tcc_numeric:.5f}s 闭式={tcc_analytic:.5f}s | "
                   f"δcc: 数值={dcc_numeric:.5f} 闭式={dcc_analytic:.5f}")


def check_scipy_cross(static, tol: float = 2e-3) -> Dict:
    """T2: 同一(有界)故障后轨迹, v1 定步长 RK4 vs scipy solve_ivp, 末态应接近.

    用故障后网络的有界振荡轨迹对比(不用发散的故障轨迹, 以免指数放大数值差).
    """

    from scipy.integrate import solve_ivp
    from bcu_3m9b.dynamics import integrate_reduced, reduced_rhs
    preset, base = static.preset, static.basevalue
    # 以预故障 SEP(相对故障后网络是被扰动但在域内的点)在故障后网络上做有界振荡.
    d0 = static.prefault.sep_delta
    w0 = np.full(preset.ngen, static.postfault.sep_omegapu * base.omega_b)
    T = 0.2
    traj = integrate_reduced(T, 1e-4, static.postfault, preset, base, d0, w0, semi_rk4=False)
    x_v1 = np.r_[traj.theta[-1], traj.omega[-1]]
    sol = solve_ivp(lambda t, x: reduced_rhs(x, static.postfault.yred, preset, base),
                    [0, T], np.r_[d0, w0], rtol=1e-11, atol=1e-13)
    x_sp = sol.y[:, -1]
    e = float(np.max(np.abs(x_v1 - x_sp)))
    return _result("轨迹: v1 RK4 vs scipy solve_ivp", e < tol, e, tol,
                   f"末态最大差={e:.2e}")


# ------------------------------- 汇总 -------------------------------

def run_all() -> List[Dict]:
    """跑全部检查, 返回结果列表."""

    from bcu_3m9b import build_static_result
    static = build_static_result()
    checks = [
        lambda: check_powerflow_residual(),
        lambda: check_sep_equilibrium(static),
        lambda: check_sep_stable(static),
        lambda: check_energy_monotonic(static),
        lambda: check_limit_cases(static),
        lambda: check_scipy_cross(static),
        lambda: check_cct_sandwich(static),
        lambda: check_smib_analytic(),
    ]
    results = []
    for fn in checks:
        try:
            results.append(fn())
        except Exception as exc:  # noqa: BLE001
            results.append(_result(getattr(fn, "__name__", "check"), False, float("nan"), 0.0,
                                   f"异常: {type(exc).__name__}: {exc}"))
    return results
