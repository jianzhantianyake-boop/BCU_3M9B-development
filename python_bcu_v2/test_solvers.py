# -*- coding: utf-8 -*-
"""P1.1 验证: SciPy 求解器统一层(与 v1 数值一致 + 更稳).

使用方法: 在 python_bcu_v2/ 目录执行  python test_solvers.py
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
from bcu_v2 import solvers


def main() -> None:
    print("=" * 66)
    print("  P1.1 · SciPy 求解器统一层验证")
    print("=" * 66)
    ok_all = True
    s = build_static_result()

    # [A] scipy SEP 与 v1 SEP 数值一致.
    from bcu_3m9b.equilibrium import solve_sep
    d_v1, _, _, _, _ = solve_sep(s.preset, s.postfault, s.basevalue, np.zeros(s.preset.ngen), 0.0)
    s.postfault.sep_delta = d_v1
    d_sp, _, perr_sp, ok_sp = solvers.solve_sep_scipy(s.preset, s.postfault, s.basevalue,
                                                      np.zeros(s.preset.ngen), 0.0)
    diff = float(np.max(np.abs(d_v1 - d_sp)))
    okA = ok_sp and diff < 1e-8
    print(f"[A] scipy SEP == v1 SEP: {'PASS' if okA else 'FAIL'}  (差={diff:.2e}, 残差={np.linalg.norm(perr_sp):.2e})")
    ok_all &= okA

    # [B] benchmark: 精度一致 + 从坏初值的稳健性.
    s.postfault.sep_delta = d_v1
    bench = solvers.benchmark_sep(s, bad_offset=1.2)
    okB = bench["match"] < 1e-8 and bench["ok_scipy_badguess"]
    print(f"[B] benchmark: {'PASS' if okB else 'FAIL'}")
    print(f"      两法解一致={bench['match']:.2e}; 坏初值收敛 v1={bench['ok_v1_badguess']} "
          f"scipy={bench['ok_scipy_badguess']}; 耗时 v1={bench['time_v1']*1e3:.1f}ms "
          f"scipy={bench['time_scipy']*1e3:.1f}ms")
    ok_all &= okB

    # [C] nlsolve 统一接口两后端一致.
    f = lambda z: np.array([z[0] ** 2 - 2.0, z[1] - z[0]])
    x1, o1, _ = solvers.nlsolve(f, np.array([1.0, 1.0]), method="scipy")
    x2, o2, _ = solvers.nlsolve(f, np.array([1.0, 1.0]), method="newton")
    okC = o1 and o2 and np.max(np.abs(x1 - x2)) < 1e-6 and abs(x1[0] - np.sqrt(2)) < 1e-6
    print(f"[C] nlsolve scipy/newton 一致: {'PASS' if okC else 'FAIL'}  (解差={np.max(np.abs(x1-x2)):.2e})")
    ok_all &= okC

    print("-" * 66)
    print(f"  总计: {'全部通过' if ok_all else '有未通过项'}")


if __name__ == "__main__":
    main()
