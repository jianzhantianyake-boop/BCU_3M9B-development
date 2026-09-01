# -*- coding: utf-8 -*-
"""P1.3 验证: 严格 DAE 级 SPM(连续法) 的一致性与稳健性.

使用方法: 在 python_bcu_v2/ 目录执行  python test_spm_dae.py
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
from bcu_3m9b import build_static_result
from bcu_3m9b.dynamics import integrate_reduced
from bcu_3m9b.spm import simulate_spm
from bcu_v2.spm_dae import simulate_spm_dae


def main() -> None:
    print("=" * 68)
    print("  P1.3 · 严格 DAE 级 SPM(连续法) 验证")
    print("=" * 68)
    ok_all = True
    s = build_static_result()
    preset, base = s.preset, s.basevalue
    sep = s.postfault.sep_delta
    w0 = np.full(preset.ngen, s.postfault.sep_omegapu * base.omega_b)

    # [A] 内部一致性: 同一 DAE 用 RK45 与 Radau 两个高精度积分器, 末端应吻合(验证 DAE 求解).
    d_init = s.prefault.sep_delta  # 相对故障后网络是被扰动但有界的点
    T = 0.3
    dae_rk = simulate_spm_dae(T, 2e-3, s.postfault, preset, base, d_init, w0, method="RK45")
    dae_rad = simulate_spm_dae(T, 2e-3, s.postfault, preset, base, d_init, w0, method="Radau")
    diff_int = float(np.max(np.abs(dae_rk["delta_coi"][-1] - dae_rad["delta_coi"][-1])))
    okA = diff_int < 1e-3
    v1 = simulate_spm(T, 2e-3, s.postfault, preset, base, d_init, w0)
    diff_v1 = float(np.max(np.abs(v1["delta_coi"][-1] - dae_rk["delta_coi"][-1])))
    print(f"[A] DAE 内部一致(RK45 vs Radau): {'PASS' if okA else 'FAIL'}  (最大差={diff_int:.2e})")
    print(f"      (参考: DAE vs v1 固定步 差={diff_v1:.2e}, 属方法/步长差异, 正常)")
    ok_all &= okA

    # [B] 稳健性: 在若干"物理上稳定但 v1 数值易失败"的清除态上比成功率.
    #     清除态由约简故障积分得到; 这些小故障(<CCT~0.24s)物理上稳定.
    d0 = s.prefault.sep_delta
    w0p = np.full(preset.ngen, s.prefault.sep_omegapu * base.omega_b)
    clear_times = [0.02, 0.025, 0.03, 0.05, 0.08, 0.10, 0.12]
    v1_ok = dae_ok = 0
    detail = []
    for tc in clear_times:
        dgen = integrate_reduced(tc, 1e-3, s.fault, preset, base, d0, w0p).theta[-1]
        # v1: 抛异常即失败
        try:
            r1 = simulate_spm(0.5, 2e-3, s.postfault, preset, base, dgen, w0)
            s1 = bool(np.all(np.isfinite(r1["delta_coi"])))
        except Exception:
            s1 = False
        # DAE: success 且有界
        try:
            r2 = simulate_spm_dae(0.5, 2e-3, s.postfault, preset, base, dgen, w0)
            s2 = bool(r2["success"] and np.all(np.isfinite(r2["delta_coi"]))
                      and np.max(np.abs(r2["delta_coi"] - sep)) < 2 * np.pi)
        except Exception:
            s2 = False
        v1_ok += int(s1); dae_ok += int(s2)
        detail.append(f"tc={tc:.3f}: v1={'ok' if s1 else 'FAIL'}, dae={'ok' if s2 else 'FAIL'}")
    okB = dae_ok >= v1_ok and dae_ok == len(clear_times)
    print(f"[B] 稳健性(小故障应全稳): {'PASS' if okB else 'FAIL'}  "
          f"(v1={v1_ok}/{len(clear_times)}, DAE={dae_ok}/{len(clear_times)})")
    for d in detail:
        print(f"      {d}")
    ok_all &= okB

    # [C] 刚性求解器 Radau 也能跑.
    try:
        rad = simulate_spm_dae(0.2, 2e-3, s.postfault, preset, base, d_init, w0, method="Radau")
        okC = bool(rad["success"])
    except Exception as exc:  # noqa: BLE001
        okC = False
        print("      Radau 异常:", exc)
    print(f"[C] 刚性 Radau 积分可用: {'PASS' if okC else 'FAIL'}")
    ok_all &= okC

    print("-" * 68)
    print(f"  总计: {'全部通过' if ok_all else '有未通过项'}")


if __name__ == "__main__":
    main()
