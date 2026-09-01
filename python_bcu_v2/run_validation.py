# -*- coding: utf-8 -*-
"""P0 + P1.2 一键验证: 跑全部不变量/金标准检查并打印报告, 附精确 CCT 对照.

使用方法:
    在 python_bcu_v2/ 目录下执行:  python run_validation.py
    需要 numpy + scipy(已装). 复用兄弟目录 python_bcu 的 v1(不改 v1).
"""

from pathlib import Path
import sys

# Windows 控制台切 UTF-8, 避免中文乱码.
if sys.platform == "win32":
    try:
        import ctypes
        ctypes.windll.kernel32.SetConsoleOutputCP(65001)
        ctypes.windll.kernel32.SetConsoleCP(65001)
    except Exception:
        pass
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8")
        except Exception:
            pass

# 让 `import bcu_v2` 可用(本文件在 python_bcu_v2/ 下, bcu_v2 是其子包).
ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from bcu_v2 import invariants, cct  # noqa: E402


def main() -> None:
    """跑验证套件并打印报告."""

    print("=" * 74)
    print("  BCU 平台 · P0 正确性验证 + P1.2 精确 CCT (zero-MATLAB)")
    print("=" * 74)

    results = invariants.run_all()
    npass = sum(r["passed"] for r in results)
    print(f"\n[T1 不变量 + T2 独立参照]  通过 {npass}/{len(results)}\n")
    print(f"  {'检查项':<34}{'结果':<6}{'误差':<12}{'容差':<10}")
    print("  " + "-" * 68)
    for r in results:
        flag = "PASS" if r["passed"] else "FAIL"
        err = "nan" if r["error"] != r["error"] else f"{r['error']:.2e}"
        print(f"  {r['name']:<34}{flag:<6}{err:<12}{r['tol']:<10.0e}")
        if r["detail"]:
            print(f"      └ {r['detail']}")

    # P1.2 精确 CCT 与 v1 网格 REA / 能量法 LEA 的直接对照.
    print("\n[P1.2 精确 CCT 对照 · 3 机网络约简]")
    from bcu_3m9b import build_static_result
    from bcu_v2.fixes import run_experiment_clean
    static = build_static_result()
    res = run_experiment_clean(static, cct_samples=21)
    precise, found = cct.precise_cct_reduced(static, tol=5e-5)
    print(f"  能量法 LEA CCT        = {res['lea'].cct:.5f} s")
    print(f"  v1 网格 REA CCT(21)   = {res['rea_cct']:.5f} s   (固定网格, 粗)")
    print(f"  P1.2 事件驱动精确 CCT = {precise:.5f} s   (二分+事件, 细; found={found})")

    print("\n[结论]")
    if npass == len(results):
        print("  全部检查通过: 平台在物理不变量与 SMIB 等面积金标准下自洽, 且精确 CCT 收敛.")
    else:
        print("  有检查未通过, 见上表 detail; 这是发现 v1 或方法问题的信号(可能是 v1 待修项).")


if __name__ == "__main__":
    main()
