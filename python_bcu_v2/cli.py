# -*- coding: utf-8 -*-
"""BCU v2 命令行入口(配置驱动) —— 镜像 MATLAB 的 run_bcu.m.

用法示例:
    python cli.py show                         # 打印当前配置摘要
    python cli.py list                         # 列出可选 mode / case
    python cli.py run                           # 用 config.yaml 跑 mode
    python cli.py run --mode reduced_region     # 临时覆盖 mode
    python cli.py run --case case39_modified --auto-params   # 切 39 母线
    python cli.py run --fault-line 9,6 --tunit 5e-4 --grid 15
    python cli.py run --set Pm=[0.9,1.3,0.95]   # 通用键值覆盖
    python cli.py validate                      # P0 物理不变量 + 金标准
    python cli.py xval                          # T3 与 MATLAB 交叉验证
"""

from pathlib import Path
import argparse
import inspect
import json
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

from bcu_v2 import config as C  # noqa: E402


def _parse_value(s: str):
    """把命令行字符串解析成 Python 值(优先 JSON, 退回原字符串)."""

    try:
        return json.loads(s)
    except Exception:
        return s


def _collect_overrides(args) -> dict:
    """从命令行参数收集配置覆盖 dict."""

    ov = {}
    if args.mode is not None:
        ov["mode"] = args.mode
    if args.case is not None:
        ov["case"] = args.case
    if args.auto_params:
        ov["auto_params"] = True
    if args.fault_line is not None:
        ov["faultline"] = [int(x) for x in args.fault_line.split(",")]
    if args.tfault is not None:
        ov["Tfault"] = float(args.tfault)
    if args.tunit is not None:
        ov["Tunit"] = float(args.tunit)
    if args.grid is not None:
        ov["region_grid"] = int(args.grid)
    for item in (args.set or []):
        if "=" not in item:
            print(f"[warn] 忽略非法 --set '{item}'(应为 key=value)"); continue
        k, v = item.split("=", 1)
        ov[k.strip()] = _parse_value(v.strip())
    return ov


def _load_cfg(args) -> dict:
    cfg = C.load_config(getattr(args, "config", None))
    return C.apply_overrides(cfg, _collect_overrides(args))


def cmd_show(args) -> None:
    C.print_summary(_load_cfg(args))


def cmd_list(args) -> None:
    print("可选 mode:")
    for m in C.MODES:
        print(f"  - {m}")
    print("可选 case:", ", ".join(C.CASES))


def cmd_validate(args) -> None:
    from bcu_v2 import invariants
    results = invariants.run_all()
    npass = sum(r["passed"] for r in results)
    print(f"P0 验证: 通过 {npass}/{len(results)}")
    for r in results:
        print(f"  [{'PASS' if r['passed'] else 'FAIL'}] {r['name']}  误差={r['error']:.2e}")


def cmd_xval(args) -> None:
    from bcu_v2 import matlab_xval
    path = matlab_xval.default_baseline_path()
    if not path.exists():
        print(f"未找到 MATLAB 参考: {path}"); return
    results = matlab_xval.run_xval(path)
    npass = sum(r["passed"] for r in results)
    print(f"T3 交叉验证: 通过 {npass}/{len(results)}")
    for r in results:
        print(f"  [{'PASS' if r['passed'] else 'FAIL'}] {r['layer']}  误差={r['error']:.2e}"
              + (f"  ({r['detail']})" if r["detail"] else ""))


def cmd_run(args) -> None:
    cfg = _load_cfg(args)

    print("==================== 配置校验 ====================")
    ok, msgs = C.validate_config(cfg)
    for m in msgs:
        print("  " + m)
    if not ok:
        print("配置校验未通过, 请修正 config.yaml 或命令行参数后重试.")
        sys.exit(1)

    C.print_summary(cfg)
    if cfg.get("save_snapshot", True):
        jp = C.save_snapshot(cfg)
        print(f"  配置快照已存: {jp}")

    from bcu_3m9b import experiments
    func = experiments.MODES[cfg["mode"]]

    # 按配置构建 static(两机模式不需要), 并按 mode 签名注入相关参数.
    static = C.build_static_from_config(cfg)
    candidate = {"grid_points": int(cfg["region_grid"]),
                 "postfault_time": float(cfg["postfault_time"]),
                 "tunit": float(cfg["Tunit"])}
    sig = inspect.signature(func)
    kwargs = {k: v for k, v in candidate.items() if k in sig.parameters}

    print(f"\n---------- 运行 mode: {cfg['mode']} ----------")
    func(static, **kwargs)

    # 运行后自检(reduced_cct/numerical): SEP 残差 + 能量法 LEA<=时域 REA(均用 v2, 任意 ngen).
    if cfg.get("run_selfcheck", True) and cfg["mode"] in ("reduced_cct", "reduced_numerical") and static is not None:
        print("\n---------- 运行后自检 ----------")
        from bcu_v2 import cct as _cct
        from bcu_v2 import cuep as _cuep
        import numpy as np
        se = float(np.linalg.norm(static.postfault.sep_perr))
        rea, _ = _cct.precise_cct_reduced(static, tol=5e-5)
        print(f"  postfault SEP 残差 : {se:.2e}  ({'PASS' if se < 1e-6 else 'FAIL'})")
        # 通用能量法 CUEP(closest-UEP, 任意 ngen); 9 与 39 母线均能出 LEA.
        lea_res = _cuep.energy_lea_cct(static)
        if lea_res.found:
            cr = lea_res.cuep
            print(f"  能量法 CUEP        : type-1, V(CUEP)={cr.v_cuep:.4f}, 候选 type-1 UEP={cr.n_type1}")
            print(f"  LEA / REA CCT      : {lea_res.lea:.4f} / {rea:.4f} s  (ngen={static.preset.ngen})")
            print(f"  能量法保守 LEA<=REA: {'PASS' if lea_res.lea <= rea + 1e-3 else 'FAIL'}")
        else:
            print(f"  时域精确 REA CCT   : {rea:.4f} s  (ngen={static.preset.ngen}; "
                  f"能量法 CUEP 未找到: {lea_res.note})")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(prog="cli.py", description="BCU v2 配置驱动命令行入口")
    sub = p.add_subparsers(dest="cmd")

    def add_common(sp):
        sp.add_argument("--config", help="配置文件路径(默认 config.yaml)")
        sp.add_argument("--mode"); sp.add_argument("--case")
        sp.add_argument("--auto-params", action="store_true", dest="auto_params")
        sp.add_argument("--fault-line"); sp.add_argument("--tfault")
        sp.add_argument("--tunit"); sp.add_argument("--grid")
        sp.add_argument("--set", action="append", help="通用覆盖 key=value(可多次)")

    for name, fn in (("show", cmd_show), ("run", cmd_run)):
        sp = sub.add_parser(name); add_common(sp); sp.set_defaults(func=fn)
    for name, fn in (("list", cmd_list), ("validate", cmd_validate), ("xval", cmd_xval)):
        sp = sub.add_parser(name); sp.set_defaults(func=fn)

    args = p.parse_args(argv)
    if not getattr(args, "cmd", None):
        p.print_help(); return
    # 为无 --config 的子命令补默认 None.
    if not hasattr(args, "config"):
        args.config = None
    args.func(args)


if __name__ == "__main__":
    main()
