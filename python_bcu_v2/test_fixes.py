# -*- coding: utf-8 -*-
"""P0.2 正确性回归: 验证 v1 的 5 个隐患**已就地修复**(不再对比旧错误行为)。

背景:
    这 5 个隐患原先只在 v2 的 fixes.py 里绕开(对比 v1 旧错误); 现已回灌到冻结解除后的 v1
    (equilibrium.solve_sep 加 inplace / bcu.run_bcu_experiment 用 inplace=False /
    dynamics.find_exitpoint 排除伪过零 / dynamics.trajectory_stable 默认有界性判据 /
    spm.solve_algebraic 用 scipy.root 稳健化)。本脚本改为验证 **v1 自身已正确**, 并与 v2 的
    fixes.* 参考实现对齐。

使用方法: 在 python_bcu_v2/ 目录执行  python test_fixes.py
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
from bcu_3m9b.bcu import run_bcu_experiment
from bcu_3m9b.dynamics import integrate_reduced, find_exitpoint, trajectory_stable
from bcu_3m9b.equilibrium import solve_sep
from bcu_3m9b.spm import solve_algebraic
from bcu_v2 import fixes


def main() -> None:
    print("=" * 66)
    print("  P0.2 · v1 隐患已回灌修复 · 正确性回归")
    print("=" * 66)
    ok_all = True

    # 隐患1(已修): v1 solve_sep 提供 inplace=False, 求异地平衡点不改写 state.sep_*.
    s = build_static_result()
    before = s.postfault.sep_delta.copy()
    solve_sep(s.preset, s.postfault, s.basevalue,
              delta0=np.zeros(s.preset.ngen) + 0.3, omega0=0.0, inplace=False)
    ok1 = np.array_equal(before, s.postfault.sep_delta)
    print(f"[1] v1 solve_sep(inplace=False) 无副作用: {'PASS' if ok1 else 'FAIL'}")
    ok_all &= ok1

    # 隐患2(已修): v1 run_bcu_experiment 求 CUEP 用 inplace=False, 不再污染 postfault.SEP.
    sa = build_static_result(); sep_a = sa.postfault.sep_delta.copy()
    run_bcu_experiment(sa)
    v1_pollute = float(np.max(np.abs(sep_a - sa.postfault.sep_delta)))
    ok2 = v1_pollute == 0.0
    print(f"[2] v1 run_bcu_experiment 不污染 postfault.SEP: {'PASS' if ok2 else 'FAIL'}  "
          f"(污染={v1_pollute:.1e})")
    ok_all &= ok2

    # 隐患3(已修): v1 find_exitpoint 现已排除初始伪过零(index>0), 与 fixes 参考一致.
    s3 = build_static_result()
    d0 = s3.prefault.sep_delta
    w0 = np.full(s3.preset.ngen, s3.prefault.sep_omegapu * s3.basevalue.omega_b)
    traj = integrate_reduced(0.2, 1e-3, s3.fault, s3.preset, s3.basevalue, d0, w0)
    idx_v1 = find_exitpoint(traj, s3.postfault, s3.preset)
    idx_fix = fixes.find_exitpoint_fixed(traj, s3.postfault, s3.preset)
    ok3 = idx_v1 > 0 and idx_v1 == idx_fix
    print(f"[3] v1 find_exitpoint 排除伪过零 & 与参考一致: {'PASS' if ok3 else 'FAIL'}  "
          f"(v1 index={idx_v1}, fixes index={idx_fix})")
    ok_all &= ok3

    # 隐患4(已修): v1 trajectory_stable 默认有界性判据, 轻阻尼下正确判稳.
    s4 = build_static_result()
    w0b = np.full(s4.preset.ngen, s4.postfault.sep_omegapu * s4.basevalue.omega_b)
    stable_traj = integrate_reduced(1.0, 1e-3, s4.postfault, s4.preset, s4.basevalue,
                                    s4.prefault.sep_delta, w0b)
    ok4_v1 = trajectory_stable(stable_traj, s4.postfault, s4.preset)  # 默认 criterion='bounded'
    ok4_ref = fixes.is_stable_bounded(stable_traj, s4.postfault.sep_delta)
    ok4 = bool(ok4_v1 and ok4_v1 == ok4_ref)
    print(f"[4] v1 trajectory_stable(有界性) 轻阻尼判稳 & 与参考一致: {'PASS' if ok4 else 'FAIL'}")
    ok_all &= ok4

    # 隐患5(已修): v1 solve_algebraic 用 scipy.root + 多初值稳健化, 冷启动成功率显著.
    # 诚实说明: 仍有个别病态清除态冷启动不收敛(单纯换求解器非万能), 真正根治靠 v2 spm_dae 的
    # 连续法(见 test_spm_dae.py: DAE 7/7)。此处断言 scipy 冷启动成功率 >= 5/6(达到原 v1 热启动水平).
    s5 = build_static_result()
    v1_ok = n = 0
    for tc in (0.05, 0.10, 0.15, 0.20, 0.25, 0.30):
        dgen = integrate_reduced(tc, 1e-3, s5.fault, s5.preset, s5.basevalue,
                                 s5.prefault.sep_delta, w0).theta[-1]
        _, o1, _ = solve_algebraic(dgen, s5.postfault, s5.preset)  # 冷启动(不传 guess)
        n += 1; v1_ok += int(o1)
    ok5 = v1_ok >= 5
    print(f"[5] v1 solve_algebraic 冷启动稳健(>=5/6, 余病态态需 DAE 连续法): "
          f"{'PASS' if ok5 else 'FAIL'}  (v1={v1_ok}/{n})")
    ok_all &= ok5

    print("-" * 66)
    print(f"  总计: {'全部通过' if ok_all else '有未通过项'}")
    return 0 if ok_all else 1


if __name__ == "__main__":
    raise SystemExit(main())
