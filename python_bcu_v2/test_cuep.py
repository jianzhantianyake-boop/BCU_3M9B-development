# -*- coding: utf-8 -*-
"""通用能量法 CUEP + LEA CCT 验证(bcu_v2.cuep, 任意 ngen 的 closest-UEP 法).

覆盖:
    [1] 9 母线 CUEP vs MATLAB baseline(交叉验证, 应 ~1e-3 内);
    [2] 9 母线 LEA vs MATLAB(0.2274 vs 0.2275);
    [3] 9 母线 通用法 CUEP == 3 机专用 reconstruct_cuep(方法内部一致);
    [4] 9 母线 CUEP 物理不变量: type-1 + V(CUEP)>V(SEP) + LEA<=REA;
    [5] 39 母线 能出 CUEP: type-1 + V(CUEP)>V(SEP)(平台由此从 3 机推广到多机);
    [6] 39 母线 能量法保守: 0 < LEA <= REA(时域真值), 即 39 母线也能出 LEA CCT.

运行: python test_cuep.py
"""

import sys
from pathlib import Path

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
for _p in (str(ROOT), str(ROOT.parent / "python_bcu")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np

from bcu_v2 import cuep as CU
from bcu_v2 import cct as CCT
from bcu_v2 import config as C
from bcu_v2 import matlab_xval as MX


def _build(case, auto=False, faultline=None):
    cfg = C.load_config()
    ov = {"case": case}
    if auto:
        ov["auto_params"] = True
    if faultline is not None:
        ov["faultline"] = faultline
    return C.build_static_from_config(C.apply_overrides(cfg, ov))


def main() -> int:
    print("=" * 62)
    print("  通用能量法 CUEP + LEA CCT 验证 (closest-UEP, 任意 ngen)")
    print("=" * 62)
    results = []

    # --------------------------- 9 母线 ---------------------------
    s9 = _build("case9_v2")
    cres9 = CU.controlling_uep(s9)
    lea9 = CU.energy_lea_cct(s9)
    rea9, _ = CCT.precise_cct_reduced(s9, tol=5e-5)

    # [1] CUEP vs MATLAB baseline.
    import scipy.io as sio
    bpath = MX.default_baseline_path()
    if bpath.exists():
        d = sio.loadmat(str(bpath), squeeze_me=True, struct_as_record=False)
        cuep_mat = np.asarray(d["postfault"].CUEP_delta, dtype=float)
        e_cuep = float(np.max(np.abs(cres9.cuep - cuep_mat)))
        ok1 = e_cuep < 5e-3
        print(f"[1] 9母线 CUEP vs MATLAB baseline: {'PASS' if ok1 else 'FAIL'}  (误差={e_cuep:.2e})")
        print(f"      py={np.array2string(cres9.cuep, precision=4)}")
        results.append(ok1)
        # [2] LEA vs MATLAB.
        lea_mat = float(d["Critical"].LEA.CCT)
        e_lea = abs(lea9.lea - lea_mat)
        ok2 = e_lea < 5e-3
        print(f"[2] 9母线 LEA vs MATLAB: {'PASS' if ok2 else 'FAIL'}  (py={lea9.lea:.4f}s, matlab={lea_mat:.4f}s)")
        results.append(ok2)
    else:
        print("[1][2] 跳过(未找到 MATLAB baseline)")

    # [3] 通用法 vs 3 机专用 reconstruct_cuep.
    cuep_3d, found3 = MX.reconstruct_cuep(s9)
    ok3 = bool(found3 and np.max(np.abs(cres9.cuep - cuep_3d)) < 5e-3)
    print(f"[3] 9母线 通用法 CUEP == 3机网格法: {'PASS' if ok3 else 'FAIL'}  "
          f"(差={np.max(np.abs(cres9.cuep - cuep_3d)):.2e})")
    results.append(ok3)

    # [4] 9 母线不变量.
    npos9 = int(np.sum(cres9.eig_reduced > 1e-6))
    ok4 = bool(cres9.found and npos9 == 1 and cres9.v_cuep > 0 and lea9.lea <= rea9 + 1e-3)
    print(f"[4] 9母线 不变量(type-1 + V(CUEP)>V(SEP) + LEA<=REA): {'PASS' if ok4 else 'FAIL'}")
    print(f"      正特征值数={npos9}, V(CUEP)={cres9.v_cuep:.4f}, LEA={lea9.lea:.4f}<=REA={rea9:.4f}")
    results.append(ok4)

    # --------------------------- 39 母线 ---------------------------
    s39 = _build("case39_modified", auto=True, faultline=[16, 17])
    cres39 = CU.controlling_uep(s39)
    ok5 = bool(cres39.found and int(np.sum(cres39.eig_reduced > 1e-6)) == 1 and cres39.v_cuep > 0)
    print(f"[5] 39母线 能出 CUEP(type-1 + V>V(SEP)): {'PASS' if ok5 else 'FAIL'}  "
          f"(ngen={s39.preset.ngen}, 候选 type-1 UEP={cres39.n_type1}, V(CUEP)={cres39.v_cuep:.4f})")
    results.append(ok5)

    lea39 = CU.energy_lea_cct(s39)
    rea39, _ = CCT.precise_cct_reduced(s39, tol=5e-5)
    ok6 = bool(lea39.found and 0.0 < lea39.lea <= rea39 + 1e-3)
    print(f"[6] 39母线 能量法保守(0 < LEA <= REA): {'PASS' if ok6 else 'FAIL'}  "
          f"(LEA={lea39.lea:.4f}s, REA={rea39:.4f}s)")
    results.append(ok6)

    print("-" * 62)
    npass = sum(results)
    print(f"  总计: {'全部通过' if npass == len(results) else f'{npass}/{len(results)} 通过'}")
    return 0 if npass == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
