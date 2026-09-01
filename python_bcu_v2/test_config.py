# -*- coding: utf-8 -*-
"""配置系统 + T3 的自动化回归测试.

使用方法: 在 python_bcu_v2/ 目录执行  python test_config.py
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
from bcu_v2 import config as C, matlab_xval


def main() -> None:
    print("=" * 62)
    print("  配置系统 + T3 回归测试")
    print("=" * 62)
    ok_all = True

    # [1] 加载 + 校验默认配置.
    cfg = C.load_config()
    ok, _ = C.validate_config(cfg)
    ok1 = ok and cfg["mode"] in C.MODES and cfg["case"] in C.CASES
    print(f"[1] 加载并校验 config.yaml: {'PASS' if ok1 else 'FAIL'}")
    ok_all &= ok1

    # [2] 命令行覆盖生效.
    cfg2 = C.apply_overrides(cfg, {"mode": "reduced_region", "Tunit": 5e-4})
    ok2 = cfg2["mode"] == "reduced_region" and cfg2["Tunit"] == 5e-4 and cfg["mode"] == "reduced_cct"
    print(f"[2] apply_overrides 生效且不改原 cfg: {'PASS' if ok2 else 'FAIL'}")
    ok_all &= ok2

    # [3] 9 母线: 按配置构建 static, SEP 收敛.
    s9 = C.build_static_from_config(cfg)
    ok3 = s9 is not None and s9.preset.ngen == 3 and np.linalg.norm(s9.postfault.sep_perr) < 1e-6
    print(f"[3] 9母线 static(SEP 残差<1e-6): {'PASS' if ok3 else 'FAIL'}  "
          f"(残差={np.linalg.norm(s9.postfault.sep_perr):.1e})")
    ok_all &= ok3

    # [4] 39 母线: auto_params + faultline 自动切换, SEP 收敛.
    cfg39 = C.apply_overrides(cfg, {"case": "case39_modified", "auto_params": True})
    s39 = C.build_static_from_config(cfg39)
    ok4 = s39 is not None and s39.preset.ngen == 10 and np.linalg.norm(s39.postfault.sep_perr) < 1e-6
    print(f"[4] 39母线10机 static(SEP 残差<1e-6): {'PASS' if ok4 else 'FAIL'}  "
          f"(残差={np.linalg.norm(s39.postfault.sep_perr):.1e})")
    ok_all &= ok4

    # [5] T3 与 MATLAB 交叉验证全通过.
    path = matlab_xval.default_baseline_path()
    if path.exists():
        res = matlab_xval.run_xval(path)
        ok5 = all(r["passed"] for r in res)
        print(f"[5] T3 MATLAB 交叉验证: {'PASS' if ok5 else 'FAIL'}  ({sum(r['passed'] for r in res)}/{len(res)})")
    else:
        ok5 = True
        print("[5] T3: 跳过(未找到 MATLAB 参考 .mat)")
    ok_all &= ok5

    print("-" * 62)
    print(f"  总计: {'全部通过' if ok_all else '有未通过项'}")


if __name__ == "__main__":
    main()
