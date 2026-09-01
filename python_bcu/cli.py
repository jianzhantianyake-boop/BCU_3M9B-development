# -*- coding: utf-8 -*-
"""BCU v1 命令行入口(配置驱动) —— 镜像 MATLAB 的 run_bcu.m / v2 的 cli.py。

用法示例:
    python cli.py show                          # 打印当前配置摘要
    python cli.py list                          # 列出可选 mode / case
    python cli.py run                           # 用 config.yaml 跑 mode(含运行后自检)
    python cli.py run --mode reduced_region --grid 15   # 临时覆盖
    python cli.py run --fault-line 8,9 --tunit 5e-4
    python cli.py run --set Pm=[0.9,1.3,0.95]   # 通用键值覆盖
    python cli.py validate                      # 潮流/SEP 残差 + 能量法 LEA<=REA

另有交互菜单入口:  python main.py                # 8 模式菜单(探索用)
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
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from bcu_3m9b import config as C  # noqa: E402


def _parse_value(s: str):
    """把命令行字符串解析成 Python 值(优先 JSON, 退回原字符串)。"""

    try:
        return json.loads(s)
    except Exception:
        return s


def _collect_overrides(args) -> dict:
    """从命令行参数收集配置覆盖 dict。"""

    ov = {}
    if args.mode is not None:
        ov["mode"] = args.mode
    if args.case is not None:
        ov["case"] = args.case
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
            print(f"[warn] ignore invalid --set '{item}' (expect key=value)"); continue
        k, v = item.split("=", 1)
        ov[k.strip()] = _parse_value(v.strip())
    return ov


def _load_cfg(args) -> dict:
    cfg = C.load_config(getattr(args, "config", None))
    return C.apply_overrides(cfg, _collect_overrides(args))


def cmd_show(args) -> None:
    C.print_summary(_load_cfg(args))


def cmd_list(args) -> None:
    print("available modes:")
    for m in C.MODES:
        print(f"  - {m}")
    print("available cases:", ", ".join(C.CASES))


def cmd_validate(args) -> None:
    """基础正确性自检: 潮流残差 + 预/故障后 SEP 残差 + 能量法 LEA<=REA。"""

    import numpy as np
    from bcu_3m9b import build_static_result
    from bcu_3m9b.cuep import energy_lea_cct
    from bcu_3m9b.experiments import robust_reduced_cct

    static = build_static_result()
    checks = []
    se_pre = float(np.linalg.norm(static.prefault.sep_perr))
    se_post = float(np.linalg.norm(static.postfault.sep_perr))
    checks.append(("pre-fault SEP residual", se_pre, se_pre < 1e-6))
    checks.append(("post-fault SEP residual", se_post, se_post < 1e-6))
    lea = energy_lea_cct(static)
    rea = robust_reduced_cct(static)
    checks.append(("energy-method LEA found", lea.lea, lea.found))
    checks.append(("conservative LEA<=REA", rea - lea.lea, lea.found and lea.lea <= rea + 1e-3))
    npass = sum(c[2] for c in checks)
    print(f"validate: {npass}/{len(checks)} passed")
    for name, val, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}  = {val:.3e}")
    if lea.found:
        print(f"  LEA = {lea.lea:.4f}s, REA = {rea:.4f}s, V(CUEP) = {lea.critical_energy:.4f}")


def cmd_run(args) -> None:
    cfg = _load_cfg(args)

    print("==================== config validation ====================")
    ok, msgs = C.validate_config(cfg)
    for m in msgs:
        print("  " + m)
    if not ok:
        print("config invalid; fix config.yaml or CLI args and retry.")
        sys.exit(1)

    C.print_summary(cfg)
    if cfg.get("save_snapshot", True):
        jp = C.save_snapshot(cfg)
        print(f"  config snapshot saved: {jp}")

    from bcu_3m9b import experiments
    func = experiments.MODES[cfg["mode"]]

    static = C.build_static_from_config(cfg)
    # 稳定域模式按网格密度注入(reduced_region 用 region_grid, spm_region 用 spm_region_grid).
    candidate = {"grid_points": int(cfg["spm_region_grid"] if cfg["mode"] == "spm_region"
                                    else cfg["region_grid"])}
    sig = inspect.signature(func)
    kwargs = {k: v for k, v in candidate.items() if k in sig.parameters}

    print(f"\n---------- run mode: {cfg['mode']} ----------")
    func(static, **kwargs)

    # 运行后自检(reduced_cct/numerical): SEP 残差 + 能量法 LEA<=REA.
    if cfg.get("run_selfcheck", True) and cfg["mode"] in ("reduced_cct", "reduced_numerical") and static is not None:
        print("\n---------- post-run self-check ----------")
        import numpy as np
        from bcu_3m9b.cuep import energy_lea_cct
        from bcu_3m9b.experiments import robust_reduced_cct
        se = float(np.linalg.norm(static.postfault.sep_perr))
        print(f"  post-fault SEP residual : {se:.2e}  ({'PASS' if se < 1e-6 else 'FAIL'})")
        lea = energy_lea_cct(static)
        rea = robust_reduced_cct(static)
        if lea.found:
            cr = lea.cuep
            print(f"  energy-method CUEP      : type-1, V(CUEP)={cr.v_cuep:.4f}, type-1 UEPs={cr.n_type1}")
            print(f"  LEA / REA CCT           : {lea.lea:.4f} / {rea:.4f} s")
            print(f"  conservative LEA<=REA   : {'PASS' if lea.lea <= rea + 1e-3 else 'FAIL'}")
        else:
            print(f"  REA CCT                 : {rea:.4f} s  (CUEP not found: {lea.note})")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(prog="cli.py", description="BCU v1 config-driven CLI")
    sub = p.add_subparsers(dest="cmd")

    def add_common(sp):
        sp.add_argument("--config", help="config file path (default config.yaml)")
        sp.add_argument("--mode"); sp.add_argument("--case")
        sp.add_argument("--fault-line"); sp.add_argument("--tfault")
        sp.add_argument("--tunit"); sp.add_argument("--grid")
        sp.add_argument("--set", action="append", help="generic override key=value (repeatable)")

    for name, fn in (("show", cmd_show), ("run", cmd_run)):
        sp = sub.add_parser(name); add_common(sp); sp.set_defaults(func=fn)
    for name, fn in (("list", cmd_list), ("validate", cmd_validate)):
        sp = sub.add_parser(name); sp.set_defaults(func=fn)

    args = p.parse_args(argv)
    if not getattr(args, "cmd", None):
        p.print_help(); return
    if not hasattr(args, "config"):
        args.config = None
    args.func(args)


if __name__ == "__main__":
    main()
