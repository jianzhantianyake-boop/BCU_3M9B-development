# -*- coding: utf-8 -*-
"""P2 验证: 通用装配层(P2.1) + ZIP 负荷(P2.3) + 单机参考模型(P2.2/P2.4).

使用方法: 在 python_bcu_v2/ 目录执行  python test_p2.py
"""

from pathlib import Path
import sys

if sys.platform == "win32":
    try:
        import ctypes
        ctypes.windll.kernel32.SetConsoleOutputCP(65001)
    except Exception:
        pass
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8")
        except Exception:
            pass

ROOT = Path(__file__).resolve().parent
for _p in (ROOT, ROOT.parent / "python_bcu"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np
from scipy.integrate import solve_ivp
from bcu_v2 import systems, loads, models


def main() -> None:
    print("=" * 70)
    print("  P2 · 建模扩展验证 (装配层 / ZIP 负荷 / 单机参考模型)")
    print("=" * 70)
    ok_all = True

    # -------- P2.1a: build_preset 复现 v1 9 母线 default_preset --------
    from bcu_3m9b.cases import case9_v2
    from bcu_3m9b.bcu import default_preset
    case9 = case9_v2()
    H9 = np.array([23.64, 6.4, 3.01]); Xd1_9 = np.array([0.0608, 0.1198, 0.1813])
    p9 = systems.build_preset(case9, H9, Xd1_9, damping=0.1)
    dp = default_preset()
    okA = (np.allclose(p9.m, dp.m, atol=1e-3) and np.allclose(p9.pmpu, dp.pmpu, atol=1e-3)
           and np.allclose(p9.epu, dp.epu) and np.allclose(p9.xd1, dp.xd1)
           and np.allclose(p9.d, dp.d, atol=1e-4))
    print(f"[P2.1a] build_preset 复现 v1 9母线参数: {'PASS' if okA else 'FAIL'}  "
          f"(m差={np.max(np.abs(p9.m-dp.m)):.1e}, pmpu差={np.max(np.abs(p9.pmpu-dp.pmpu)):.1e})")
    ok_all &= okA

    # -------- P2.1b: 39 母线动态案例跑通静态初始化(潮流->约简->SEP) --------
    try:
        case39, preset39 = systems.case39_dynamic()
        st39 = systems.build_static_dynamic(case39, preset39)
        okB = (st39.prefault.yred.shape == (10, 10)
               and np.linalg.norm(st39.postfault.sep_perr) < 1e-6)
        detail = (f"Yred={st39.prefault.yred.shape}, "
                  f"post SEP 残差={np.linalg.norm(st39.postfault.sep_perr):.2e}")
    except Exception as exc:  # noqa: BLE001
        okB = False; detail = f"异常: {type(exc).__name__}: {exc}"
    print(f"[P2.1b] 39母线10机静态初始化跑通: {'PASS' if okB else 'FAIL'}  ({detail})")
    ok_all &= okB

    # -------- P2.3: ZIP 负荷极限 --------
    P, _ = loads.zip_load_power(1.0, 1.0, 0.0, zp=(0, 0, 1))      # 恒功率 @V=1 -> P0
    Pz, _ = loads.zip_load_power(0.9, 1.0, 0.0, zp=(1, 0, 0))     # 恒阻抗 @V=0.9 -> 0.81 P0
    Pi, _ = loads.zip_load_power(0.9, 1.0, 0.0, zp=(0, 1, 0))     # 恒电流 @V=0.9 -> 0.9 P0
    okC = abs(P - 1.0) < 1e-12 and abs(Pz - 0.81) < 1e-12 and abs(Pi - 0.9) < 1e-12
    print(f"[P2.3] ZIP 负荷极限(P/Z/I): {'PASS' if okC else 'FAIL'}  "
          f"(P={P:.3f}, Z@0.9={Pz:.3f}, I@0.9={Pi:.3f})")
    ok_all &= okC

    # -------- P2.3b: ZIP(aP=1) 残差 == v1 恒功率残差 --------
    from bcu_3m9b import build_static_result
    from bcu_3m9b.spm import algebraic_residual, _load_power
    s = build_static_result()
    post = s.postfault; n = s.preset.ngen
    yorg = np.asarray(post.metadata["yorg_mod"], dtype=complex)
    transform = np.asarray(post.metadata["transform"], dtype=int)
    load_pq = _load_power(s.preset, transform[n:])
    nload = yorg.shape[0] - n
    z = np.r_[np.zeros(nload), np.ones(nload)]
    dgen = post.sep_delta
    r_v1 = algebraic_residual(z, dgen, yorg, load_pq, n)
    r_zip = loads.zip_algebraic_residual(z, dgen, yorg, load_pq, n, zp=(0, 0, 1), zq=(0, 0, 1))
    okD = np.allclose(r_v1, r_zip)
    print(f"[P2.3b] ZIP(aP=1) 残差 == v1 恒功率: {'PASS' if okD else 'FAIL'}  "
          f"(差={np.max(np.abs(r_v1-r_zip)):.1e})")
    ok_all &= okD

    # -------- P2.2: one-axis 在 Xd=Xd' 极限退化为经典 SMIB --------
    Xd1 = 0.3; Xe = 0.2; Vinf = 1.0; Efd = 1.2; Pm = 0.8; H = 3.5
    oa = models.OneAxisSMIB(H=H, D=0.0, Xd=Xd1, Xd1=Xd1, Xq=Xd1, Xe=Xe,
                            Vinf=Vinf, Pm=Pm, Efd=Efd, Tdo1=5.0)
    ws = models.OMEGA_S
    Pmax = Efd * Vinf / (Xd1 + Xe)
    T = 1.0
    x0 = np.array([0.3, ws, Efd])  # 从非平衡角出发, E'q 初值=Efd
    sol_oa = solve_ivp(lambda t, x: models.one_axis_rhs(x, oa), [0, T], x0,
                       rtol=1e-10, atol=1e-12, dense_output=True)
    # 经典 SMIB: [δ, ν=dδ/dt], M=2H/ωs, Pe=Pmax sinδ
    M = 2.0 * H / ws
    def classic(t, y):
        return [y[1], (Pm - Pmax * np.sin(y[0])) / M]
    sol_cl = solve_ivp(classic, [0, T], [0.3, 0.0], rtol=1e-10, atol=1e-12, dense_output=True)
    tt = np.linspace(0, T, 200)
    d_oa = sol_oa.sol(tt)[0]; d_cl = sol_cl.sol(tt)[0]
    eqvar = float(np.max(np.abs(sol_oa.sol(tt)[2] - Efd)))
    okE = np.max(np.abs(d_oa - d_cl)) < 1e-6 and eqvar < 1e-9
    print(f"[P2.2] one-axis 在 Xd=Xd' 退化为经典: {'PASS' if okE else 'FAIL'}  "
          f"(δ最大差={np.max(np.abs(d_oa-d_cl)):.1e}, E'q偏离Efd={eqvar:.1e})")
    ok_all &= okE

    # -------- P2.4: GFM 下垂稳态 P=Pset, 且 rhs 在 δss 为 0 --------
    gfm = models.GFMDroopSMIB(E=1.0, Vinf=1.0, X=0.5, Pset=0.8, mp=0.05)
    dss = models.gfm_equilibrium(gfm)
    Pss = models.gfm_power(dss, gfm)
    rhs0 = float(models.gfm_rhs([dss], gfm)[0])
    okF = abs(Pss - gfm.Pset) < 1e-12 and abs(rhs0) < 1e-12
    print(f"[P2.4] GFM 下垂稳态 P=Pset: {'PASS' if okF else 'FAIL'}  "
          f"(δss={dss:.4f}, P(δss)={Pss:.4f}, rhs={rhs0:.1e})")
    ok_all &= okF

    print("-" * 70)
    print(f"  总计: {'全部通过' if ok_all else '有未通过项'}")


if __name__ == "__main__":
    main()
