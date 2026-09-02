"""运行固定种子的 SPM 联合平衡根审计并保存 JSON 报告。"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT.parent / "python_bcu") not in sys.path:
    sys.path.insert(0, str(ROOT.parent / "python_bcu"))

from bcu_3m9b import build_static_result
from bcu_v2.spm_joint_audit import audit_spm_joint_roots


def load_seed_file(path: Path) -> tuple[list[list[float]], list[dict[str, Any]]]:
    """Extract replayable MGP angle and complete network seeds from a report.

    ``spm_mgp_iterations.json`` contributes every finite ``update_delta`` and
    ``end_delta`` in its trace.  A physical network seed is accepted only when
    the report explicitly contains corrected network angles and voltages (as
    in ``spm_mgp_diagnostic.json``); stored mixed-frame angles are intentionally
    ignored.  The loader is read-only and never invents missing network values.
    """

    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"seed file must contain a JSON object: {path}")
    deltas: list[list[float]] = []
    trace = payload.get("trace", [])
    if isinstance(trace, list):
        for item in trace:
            if not isinstance(item, dict):
                continue
            for key in ("update_delta", "end_delta"):
                value = item.get(key)
                if isinstance(value, list) and value:
                    try:
                        numbers = [float(x) for x in value]
                    except (TypeError, ValueError):
                        continue
                    if all(map(math.isfinite, numbers)):
                        deltas.append(numbers)
    top_delta = payload.get("mgp_delta")
    if isinstance(top_delta, list) and top_delta:
        try:
            numbers = [float(x) for x in top_delta]
        except (TypeError, ValueError):
            numbers = []
        if numbers and all(map(math.isfinite, numbers)):
            deltas.append(numbers)

    states: list[dict[str, Any]] = []
    physical = payload.get("physical_cuep")
    if isinstance(physical, dict):
        delta = physical.get("cuep_delta")
        theta = physical.get("corrected_net_theta")
        voltage = physical.get("cuep_net_voltage")
        if all(isinstance(value, list) and value for value in (delta, theta, voltage)):
            states.append({
                "delta_gen": delta,
                "theta_net": theta,
                "voltage_net": voltage,
            })
    extra_states = payload.get("states", [])
    if isinstance(extra_states, list):
        for state in extra_states:
            if isinstance(state, dict):
                states.append(state)
    return deltas, states


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--starts", type=int, default=152)
    parser.add_argument("--seed", type=int, default=20260902)
    parser.add_argument(
        "--seed-file",
        type=Path,
        action="append",
        default=[],
        help="只读读取 MGP 角度/物理网络种子报告；可重复指定",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT.parent / "validation" / "reports" / "spm_joint_audit_latest.json",
    )
    args = parser.parse_args(argv)
    if args.starts <= 0:
        parser.error("--starts must be positive")
    seed_deltas: list[list[float]] = []
    seed_states: list[dict[str, Any]] = []
    for path in args.seed_file:
        if not path.is_file():
            parser.error(f"seed file not found: {path}")
        try:
            deltas, states = load_seed_file(path)
        except (OSError, json.JSONDecodeError, ValueError) as exc:
            parser.error(str(exc))
        seed_deltas.extend(deltas)
        seed_states.extend(states)
    report = audit_spm_joint_roots(
        build_static_result(), max_starts=args.starts, random_seed=args.seed,
        seed_deltas=seed_deltas or None,
        seed_states=seed_states or None,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "starts_evaluated": report["starts_evaluated"],
        "converged_start_count": report["converged_start_count"],
        "unique_converged_root_count": report["unique_converged_root_count"],
        "used_external_ecritical": report["used_external_ecritical"],
        "explicit_delta_seeds": len(seed_deltas),
        "explicit_state_seeds": len(seed_states),
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
