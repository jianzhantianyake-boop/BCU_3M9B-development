"""第一批 24 案例的阻尼敏感性研究驱动。

八路径门禁未全部为 ``MATLAB_XVAL_FULL`` 时只输出 BLOCKED 诊断，不运行案例、不写
零值结果。门禁通过后再按固定顺序生成 CSV、JSON manifest 和环境记录。
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path


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
                    fault_bus=bus, fault_position="from",
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
    # 实际求解阶段必须由具体验收后的单案例函数填充；此处先固定 schema 和顺序。
    payload = {"schema_version": "1.0", "status": "READY", "created_at": stamp,
               "case_count": len(cases), "cases": [asdict(case) for case in cases]}
    manifest.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["case_id", "source_revision", "reference_version", "fault_bus",
                         "fault_position", "gamma", "damping_allocation", "convergence_status",
                         "failure_reason"])
    print(json.dumps({"status": "READY", "case_count": len(cases),
                      "manifest": str(manifest), "csv": str(csv_path)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
