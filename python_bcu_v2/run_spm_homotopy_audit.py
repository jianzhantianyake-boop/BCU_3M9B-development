"""运行可重复的 SPM 联合同伦诊断。"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT.parent / "python_bcu") not in sys.path:
    sys.path.insert(0, str(ROOT.parent / "python_bcu"))

from bcu_3m9b import build_static_result
from bcu_v2.spm_homotopy import continue_spm_joint_homotopy


def _load_state(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"MGP report must contain a JSON object: {path}")
    if all(key in payload for key in ("mgp_delta", "mgp_theta", "mgp_voltage")):
        return {
            "delta_gen": payload["mgp_delta"],
            "theta_net": payload["mgp_theta"],
            "voltage_net": payload["mgp_voltage"],
        }
    physical = payload.get("physical_cuep")
    if isinstance(physical, dict):
        required = ("cuep_delta", "corrected_net_theta", "cuep_net_voltage")
        if all(key in physical for key in required):
            return {
                "delta_gen": physical["cuep_delta"],
                "theta_net": physical["corrected_net_theta"],
                "voltage_net": physical["cuep_net_voltage"],
            }
    raise ValueError(
        "MGP report does not contain a complete mgp_delta/mgp_theta/mgp_voltage state"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--seed-report", type=Path,
        default=ROOT.parent / "validation" / "reports" / "spm_mgp_iterations_with_network.json",
    )
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT.parent / "validation" / "reports" / "spm_homotopy_audit_latest.json",
    )
    args = parser.parse_args(argv)
    if args.steps <= 0:
        parser.error("--steps must be positive")
    if not args.seed_report.is_file():
        parser.error(f"seed report not found: {args.seed_report}")
    try:
        state = _load_state(args.seed_report)
        result = continue_spm_joint_homotopy(build_static_result(), state, steps=args.steps)
    except (OSError, json.JSONDecodeError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    payload = result.as_dict()
    payload["schema_version"] = "1.0"
    payload["seed_report_name"] = args.seed_report.name
    payload["used_external_ecritical"] = False
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "converged": result.converged,
        "lambda_reached": result.lambda_reached,
        "equilibrium_type": result.equilibrium_type,
        "e_critical": result.e_critical,
        "network_residual": result.network_residual,
        "equilibrium_residual": result.equilibrium_residual,
        "used_external_ecritical": False,
    }, ensure_ascii=False))
    return 0 if result.converged else 1


if __name__ == "__main__":
    raise SystemExit(main())
