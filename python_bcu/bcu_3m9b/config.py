# -*- coding: utf-8 -*-
"""集中配置系统: 加载 config.yaml、校验、按配置构建 static、打印摘要、存快照。

镜像 MATLAB 平台的 bcu_config.m / run_bcu.m(以及 v2 的同名模块): 单一配置 + 命令行入口, 只调参
不改核心方程。v1 只含 3 机 9 母线经典模型(case9_v2); 更大系统/39 母线在 v2。
"""

from __future__ import annotations

import datetime as _dt
import json
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

MODES = ["reduced_cct", "reduced_numerical", "reduced_region",
         "spm_cct", "spm_numerical", "spm_region",
         "two_machine_region_3d", "two_machine_region_3d_gfl"]
CASES = ["case9_v2"]

DEFAULT_CONFIG: Dict = {
    "mode": "reduced_cct", "case": "case9_v2", "f_base": 60,
    "m": [0.1254, 0.0340, 0.0160], "damping_ratio": [0.1, 0.1, 0.1],
    "Pm": [0.8980, 1.3432, 0.9419], "xd1": [0.0608, 0.1198, 0.1813],
    "E": [1.1083, 1.1071, 1.0606],
    "faultline": [9, 6], "faultposition": 0,
    "path_energy_cal": 0, "Tfault": 0.3, "Tunit": 0.001,
    "postfault_time": 2.0, "cct_samples": 21,
    "region_grid": 21, "spm_region_grid": 15,
    "run_selfcheck": True, "save_snapshot": True,
}


def _read_simple_yaml(path: Path) -> Dict:
    """极简 YAML 读取(扁平 key: value, 用 json 逐值解析), 免 PyYAML 依赖。

    仅支持标量(数/字符串/true|false)与行内列表 [..]; 去行内 # 注释。复杂 YAML 请装 PyYAML。
    """

    out: Dict = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if ":" not in line:
            continue
        key, val = line.split(":", 1)
        key, val = key.strip(), val.strip()
        if not key or not val:
            continue
        try:
            out[key] = json.loads(val)
        except Exception:
            out[key] = val.strip('"').strip("'")
    return out


def load_config(path=None) -> Dict:
    """加载配置: 从 config.yaml 读入并覆盖默认(优先 PyYAML, 缺则用内置极简解析)。

    使用方法: path 缺省用 python_bcu/config.yaml; 返回合并后的 dict。
    """

    cfg = dict(DEFAULT_CONFIG)
    p = Path(path) if path else (Path(__file__).resolve().parents[1] / "config.yaml")
    if p.exists():
        try:
            try:
                import yaml
                loaded = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
            except ImportError:
                loaded = _read_simple_yaml(p)
            cfg.update({k: v for k, v in loaded.items() if v is not None})
        except Exception as exc:  # noqa: BLE001
            print(f"[config] failed to read {p.name} ({exc}); using built-in defaults.")
    return cfg


def apply_overrides(cfg: Dict, overrides: Dict) -> Dict:
    """把命令行覆盖(dict)合并进配置, 返回新 dict。"""

    out = dict(cfg)
    out.update({k: v for k, v in overrides.items() if v is not None})
    return out


def validate_config(cfg: Dict) -> Tuple[bool, List[str]]:
    """校验配置, 返回 (是否通过, 消息列表)。"""

    msgs, ok = [], True
    if cfg["mode"] not in MODES:
        ok = False; msgs.append(f"[ERR] mode='{cfg['mode']}' invalid, choose: {MODES}")
    else:
        msgs.append(f"[OK ] mode = {cfg['mode']}")
    if cfg["case"] not in CASES:
        ok = False; msgs.append(f"[ERR] case='{cfg['case']}' invalid, choose: {CASES}")
    else:
        msgs.append(f"[OK ] case = {cfg['case']}")
    lens = {len(cfg[k]) for k in ("m", "damping_ratio", "Pm", "xd1", "E")}
    if len(lens) != 1:
        ok = False; msgs.append(f"[ERR] generator parameter vectors length mismatch: {lens}")
    elif lens != {3}:
        ok = False; msgs.append(f"[ERR] v1 only supports 3-machine case9_v2 (got len={lens.pop()})")
    else:
        msgs.append("[OK ] generator params length = 3")
    if len(cfg["faultline"]) != 2:
        ok = False; msgs.append("[ERR] faultline must be [FromBus, ToBus]")
    if float(cfg["Tunit"]) <= 0:
        ok = False; msgs.append("[ERR] Tunit must be > 0")
    if float(cfg["Tfault"]) <= 0:
        ok = False; msgs.append("[ERR] Tfault must be > 0")
    return ok, msgs


def build_static_from_config(cfg: Dict):
    """按配置构建 v1 StaticResult(潮流->约简->SEP), 仅 case9_v2。

    使用方法: 传入配置 dict, 返回 StaticResult; two_machine 模式不需要它(返回 None)。
    """

    from .bcu import build_static_result
    from .types import Preset

    if cfg["mode"].startswith("two_machine"):
        return None

    m = np.asarray(cfg["m"], dtype=float)
    d = m * np.asarray(cfg["damping_ratio"], dtype=float)
    preset = Preset(m=m, d=d, pmpu=np.asarray(cfg["Pm"], dtype=float),
                    xd1=np.asarray(cfg["xd1"], dtype=float),
                    epu=np.asarray(cfg["E"], dtype=float),
                    path_energy_cal=int(cfg["path_energy_cal"]),
                    fault_line=np.asarray(cfg["faultline"], dtype=int),
                    fault_position=int(cfg["faultposition"]))
    static = build_static_result(preset=preset)
    static.basevalue.omega_b = 2.0 * np.pi * float(cfg["f_base"])
    return static


def print_summary(cfg: Dict) -> None:
    """打印配置摘要(镜像 run_bcu 的摘要段)。"""

    print("=" * 56)
    print("  Experiment configuration")
    print("=" * 56)
    print(f"  mode              : {cfg['mode']}")
    print(f"  case / f_base     : {cfg['case']} @ {cfg['f_base']} Hz")
    print(f"  generators        : {len(cfg['m'])}")
    print(f"  inertia m         : {cfg['m']}")
    print(f"  damping d/m       : {cfg['damping_ratio']}")
    print(f"  mech power Pm     : {cfg['Pm']}")
    print(f"  fault line / pos  : {cfg['faultline']} / pos {cfg['faultposition']}")
    print(f"  PathEnergyCal/grid: {cfg['path_energy_cal']} / region {cfg['region_grid']}")
    print(f"  Tfault / Tunit    : {cfg['Tfault']} s / {cfg['Tunit']:g} s")
    print("=" * 56)


def save_snapshot(cfg: Dict, resdir=None) -> Path:
    """把配置存快照到 results/(json), 返回 json 路径。"""

    resdir = Path(resdir or (Path(__file__).resolve().parents[1] / "results"))
    resdir.mkdir(exist_ok=True)
    stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    jpath = resdir / f"snapshot_{stamp}.json"
    with open(jpath, "w", encoding="utf-8") as f:
        json.dump(cfg, f, ensure_ascii=False, indent=2)
    return jpath
