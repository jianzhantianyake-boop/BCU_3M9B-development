"""第一批 24 案例的阻尼敏感性研究驱动。

八路径门禁未全部为 ``MATLAB_XVAL_FULL`` 时只输出 BLOCKED 诊断，不运行案例、不写
零值结果。门禁通过后再按固定顺序生成 CSV、JSON manifest 和环境记录。
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import time
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "validation" / "reports" / "full_xval_latest.json"
OUT = ROOT / "experiments" / "results"


@dataclass(frozen=True)
class Case:
    case_id: str
    source_revision: str
    reference_version: str
    fault_bus: int
    fault_position: str
    gamma: float
    damping_allocation: str


RESULT_FIELDS = [
    "case_id", "source_revision", "reference_version", "fault_bus", "fault_position",
    "gamma", "damping_allocation", "reduced_e_critical", "spm_e_critical",
    "reduced_lea_cct", "spm_lea_cct", "reduced_rea_cct", "spm_rea_cct",
    "lea_rea_error_ms", "lea_rea_error_pct", "cross_model_gap_ms", "cuep_branch_id",
    "equilibrium_residual", "network_residual", "flag_cct", "convergence_status",
    "failure_reason", "runtime_seconds",
]


def _source_revision() -> str:
    """Return the integrated repository HEAD without changing its state."""

    try:
        return subprocess.check_output(
            ["git", "-c", f"safe.directory={ROOT}", "-C", str(ROOT), "rev-parse", "HEAD"],
            text=True, stderr=subprocess.STDOUT,
        ).strip()
    except Exception:
        return "unknown"


def _allocation_ratios(gamma: float, allocation: str) -> list[float]:
    """Map the registered study labels to explicit damping ratios.

    ``U`` scales all three baseline ratios; ``G1``--``G3`` scale only the
    named generator.  The mapping is recorded here so a later batch cannot
    silently reinterpret the allocation labels.
    """

    ratios = [0.1, 0.1, 0.1]
    if allocation == "U":
        return [gamma] * 3
    index = {"G1": 0, "G2": 1, "G3": 2}.get(allocation)
    if index is None:
        raise ValueError(f"unknown damping allocation: {allocation}")
    ratios[index] = gamma
    return ratios


def _run_case(case: Case) -> dict:
    """Run one case after the global gate has passed; failures remain explicit."""

    started = time.perf_counter()
    row = {field: None for field in RESULT_FIELDS}
    row.update(asdict(case), source_revision=_source_revision(),
               convergence_status="FAILED")

    def finite_or_none(value):
        try:
            return float(value) if np.isfinite(value) else None
        except (TypeError, ValueError):
            return None
    try:
        # Imports are delayed so a BLOCKED gate never initializes a solver or
        # creates runtime caches/results.
        sys.path.insert(0, str(ROOT / "python_bcu_v2"))
        sys.path.insert(0, str(ROOT / "python_bcu"))
        from bcu_v2 import config
        from bcu_v2.cuep import controlling_uep
        from bcu_v2.fixes import run_experiment_clean
        from bcu_v2.spm_cuep import spm_self_contained_cct

        position = 0 if case.fault_position == "from" else 1
        cfg = config.apply_overrides(config.load_config(), {
            "mode": "reduced_cct", "faultposition": position,
            "damping_ratio": _allocation_ratios(case.gamma, case.damping_allocation),
        })
        static = config.build_static_from_config(cfg)
        reduced = run_experiment_clean(static, cct_samples=21)
        cuep = controlling_uep(static)
        spm = spm_self_contained_cct(static)
        row["reduced_e_critical"] = float(cuep.v_cuep) if cuep.found else None
        row["reduced_lea_cct"] = float(reduced["lea"].cct)
        row["reduced_rea_cct"] = float(reduced["rea_cct"])
        row["spm_e_critical"] = finite_or_none(spm.cuep.e_critical)
        row["spm_lea_cct"] = finite_or_none(spm.cct) if spm.converged else None
        row["equilibrium_residual"] = finite_or_none(spm.cuep.equilibrium_residual)
        row["network_residual"] = finite_or_none(spm.cuep.network_residual)
        row["cuep_branch_id"] = getattr(spm.cuep, "equilibrium_type", None)
        row["flag_cct"] = bool(spm.converged)
        row["convergence_status"] = "PASSED" if spm.converged else "UNVERIFIED"
        if row["reduced_lea_cct"] is not None and row["reduced_rea_cct"] is not None:
            row["lea_rea_error_ms"] = 1000.0 * (row["reduced_lea_cct"] - row["reduced_rea_cct"])
            row["lea_rea_error_pct"] = 100.0 * (row["reduced_lea_cct"] - row["reduced_rea_cct"]) / max(abs(row["reduced_rea_cct"]), 1e-15)
        if row["spm_lea_cct"] is not None and row["reduced_lea_cct"] is not None:
            row["cross_model_gap_ms"] = 1000.0 * (row["spm_lea_cct"] - row["reduced_lea_cct"])
        if not spm.converged:
            row["failure_reason"] = spm.exit_reason or "SPM self-contained CCT unavailable"
    except Exception as exc:  # noqa: BLE001
        row["failure_reason"] = f"{type(exc).__name__}: {exc}"
    row["runtime_seconds"] = round(time.perf_counter() - started, 6)
    return row


def build_cases() -> list[Case]:
    cases = []
    for fault in ("F9", "F6"):
        bus = int(fault[1:])
        for gamma in (0.5, 1.0, 1.5):
            for allocation in ("U", "G1", "G2", "G3"):
                cases.append(Case(
                    case_id=f"{fault}-gamma{gamma:g}-{allocation}",
                    source_revision="4bc054c81b7193318ac2bf179a248dca3cc3341b",
                    reference_version="pending-matlab-compact-v1",
                    fault_bus=bus, fault_position=("from" if fault == "F9" else "to"),
                    gamma=gamma, damping_allocation=allocation))
    return cases


def gate_status() -> tuple[bool, list[str]]:
    if not REPORT.exists():
        return False, ["缺少 validation/reports/full_xval_latest.json"]
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    entries = report.get("entries", [])
    if len(entries) != 8:
        return False, [f"八路径条目数为 {len(entries)}，不是 8"]
    incomplete = [item["name"] for item in entries if item.get("status") != "MATLAB_XVAL_FULL"]
    if incomplete:
        return False, ["以下路径未完成 MATLAB_XVAL_FULL: " + ", ".join(incomplete)]
    return True, []


def main() -> int:
    ok, reasons = gate_status()
    if not ok:
        result = {
            "schema_version": "1.0",
            "status": "BLOCKED",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "case_count": 24,
            "cases": [],
            "failure_reason": reasons,
        }
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 2

    cases = build_cases()
    OUT.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%SZ")
    manifest = OUT / f"damping_{stamp}.json"
    csv_path = OUT / f"damping_{stamp}.csv"
    rows = [_run_case(case) for case in cases]
    payload = {"schema_version": "1.0", "status": "COMPLETED", "created_at": stamp,
               "case_count": len(cases), "cases": rows}
    manifest.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"status": "COMPLETED", "case_count": len(cases),
                      "manifest": str(manifest), "csv": str(csv_path)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
