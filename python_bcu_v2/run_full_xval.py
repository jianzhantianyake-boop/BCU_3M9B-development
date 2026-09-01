# -*- coding: utf-8 -*-
"""全量交叉验证: Python 各路径 vs MATLAB 平台参考(逐路径)。

覆盖(逐步补齐, 见 验证覆盖矩阵_CN.md):
    [reduced_cct]      -> run_matlab_xval.py (T3, 已有 8/8)
    [reduced_region]   平衡点集合(SEP + type-1 UEP)
    [two_machine_3d]   平衡点集合(含 D2=0.5 对齐)
    [two_machine_gfl]  平衡点集合
    [reduced_numerical] 三段轨迹末端       (待补)
    [spm_cct]          SPM 退出点/MGP/CUEP/LEA (待补, 阶段 C)

参考文件由 matlab_platform/verify/export_*.m 用 `matlab -batch` 生成(不改 B3_MM)。
运行: python run_full_xval.py
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
import scipy.io as sio

VERIFY = ROOT.parent / "matlab_platform" / "verify"


def _match_eps(xeps_m, flags_m, eps_py, pyxep, pyflag, ndim=None):
    """把 MATLAB 平衡点集合与 Python 集合逐点最近匹配, 返回 (最大坐标误差, flag 是否全一致)。"""

    max_err = 0.0
    flag_ok = True
    for xm, fm in zip(xeps_m, flags_m):
        xm = np.asarray(xm, dtype=float).reshape(-1)
        d = ndim if ndim else xm.size
        best, berr = None, np.inf
        for e in eps_py:
            xp = np.asarray(pyxep(e), dtype=float).reshape(-1)[:d]
            err = np.linalg.norm(xp - xm[:d])
            if err < berr:
                berr, best = err, e
        max_err = max(max_err, berr)
        if best is None or int(pyflag(best)) != int(fm):
            flag_ok = False
    return max_err, flag_ok


def verify_reduced_region() -> dict:
    from bcu_v2 import config as C
    from bcu_3m9b.experiments import find_reduced_equilibria
    f = VERIFY / "baseline_region.mat"
    if not f.exists():
        return {"layer": "reduced_region", "error": np.nan, "tol": 5e-3, "passed": False,
                "detail": "缺 baseline_region.mat(先跑 export_region.m)"}
    d = sio.loadmat(str(f), squeeze_me=True)
    xeps = np.atleast_2d(d["xeps"]); flags = np.atleast_1d(d["flags"])
    s = C.build_static_from_config(C.load_config())
    eps = find_reduced_equilibria(s.postfault, s.preset, grid_points=21)
    err, fok = _match_eps(xeps, flags, eps, lambda e: e["xep"], lambda e: e["flag"])
    return {"layer": "reduced_region", "error": err, "tol": 5e-3,
            "passed": bool(err < 5e-3 and fok),
            "detail": f"{len(eps)} EP, flag一致={fok}"}


def verify_two_machine(gfl: bool) -> dict:
    from bcu_3m9b.experiments import _params_two_machine, find_two_machine_equilibria_3d
    tag = "two_machine_gfl" if gfl else "two_machine_3d"
    f = VERIFY / ("baseline_twomachine_gfl.mat" if gfl else "baseline_twomachine.mat")
    if not f.exists():
        return {"layer": tag, "error": np.nan, "tol": 5e-3, "passed": False,
                "detail": "缺参考(先跑 export_twomachine.m)"}
    d = sio.loadmat(str(f), squeeze_me=True)
    raw = d["xeps"]; flags = np.atleast_1d(d["flags"])
    xeps = [np.asarray(x).reshape(-1) for x in (raw if raw.dtype == object else [raw])]
    p = _params_two_machine(gfl=gfl)
    eps = find_two_machine_equilibria_3d(p, grid_points=11)
    err, fok = _match_eps(xeps, flags, eps, lambda e: e["xep"], lambda e: e["flag"], ndim=3)
    return {"layer": tag, "error": err, "tol": 5e-3,
            "passed": bool(err < 5e-3 and fok),
            "detail": f"{len(eps)} EP, flag一致={fok}"}


def verify_reduced_numerical() -> dict:
    from bcu_v2 import config as C
    from bcu_3m9b.dynamics import integrate_reduced
    f = VERIFY / "baseline_numerical.mat"
    if not f.exists():
        return {"layer": "reduced_numerical", "error": np.nan, "tol": 1e-5, "passed": False,
                "detail": "缺 baseline_numerical.mat(先跑 export_numerical.m)"}
    d = sio.loadmat(str(f), squeeze_me=True)
    # 用 COI 相对角 thetac 对比(两边同坐标; theta 因绝对角含 COI 漂移不可直接比).
    thetac_end_m = np.atleast_1d(d["thetac_end"]).astype(float)
    s = C.build_static_from_config(C.load_config())
    d0 = np.asarray(s.prefault.sep_delta, dtype=float)
    w0 = np.full(s.preset.ngen, s.prefault.sep_omegapu * s.basevalue.omega_b)
    traj = integrate_reduced(0.24, 1e-4, s.fault, s.preset, s.basevalue, d0, w0)
    err = float(np.max(np.abs(traj.thetac[-1] - thetac_end_m)))
    return {"layer": "reduced_numerical", "error": err, "tol": 1e-5,
            "passed": bool(err < 1e-5), "detail": "故障段末端切除态 thetac(0.24s)"}


def verify_spm_cct() -> dict:
    """SPM 能量法 CCT: 势能/网络解/发电机功率/CCT 机制已交叉验证; E_crit 暂用 MATLAB 参考。

    诚实说明: 自足求 SPM CUEP 网络态的分支选择尚未闭环(见 spm_energy.py 头部与验证覆盖矩阵),
    故此处 E_critical 取 MATLAB baseline(3.3757), 验证的是"给定 E_crit, Python 独立复现 CCT"。
    """
    from bcu_v2 import config as C
    from bcu_v2 import spm_energy as SE
    f = VERIFY / "baseline_spm.mat"
    if not f.exists():
        return {"layer": "spm_cct*", "error": np.nan, "tol": 5e-3, "passed": False,
                "detail": "缺 baseline_spm.mat(先跑 export_spm.m)"}
    d = sio.loadmat(str(f), squeeze_me=True, struct_as_record=False)["ref"]
    ecrit = float(d.E_critical); cct_m = float(d.CCT_LEA)
    s = C.build_static_from_config(C.apply_overrides(C.load_config(), {"mode": "spm_cct"}))
    cct, ok = SE.spm_fault_energy_cct(s, ecrit)
    err = abs(cct - cct_m)
    return {"layer": "spm_cct*", "error": err, "tol": 5e-3, "passed": bool(ok and err < 5e-3),
            "detail": f"CCT {cct:.4f} vs {cct_m:.4f}; E_crit=3.3757 取自MATLAB(自足CUEP待做)"}


def main() -> int:
    print("=" * 66)
    print("  全量交叉验证 · Python vs MATLAB 平台(逐路径)")
    print("=" * 66)
    results = [
        verify_reduced_region(),
        verify_reduced_numerical(),
        verify_two_machine(gfl=False),
        verify_two_machine(gfl=True),
        verify_spm_cct(),
    ]
    npass = sum(r["passed"] for r in results)
    for r in results:
        tag = "PASS" if r["passed"] else "FAIL"
        print(f"  [{tag}] {r['layer']:<22} 误差={r['error']:.2e}  容差={r['tol']:.0e}  ({r['detail']})")
    print("-" * 66)
    print(f"  小计: {npass}/{len(results)} 通过")
    print("  注: reduced_cct 见 run_matlab_xval.py(T3 8/8)。spm_cct* 的 E_crit 取自 MATLAB 参考,")
    print("      自足求 SPM CUEP 网络态(分支选择)是明确的下一里程碑, 见 验证覆盖矩阵_CN.md。")
    return 0 if npass == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
