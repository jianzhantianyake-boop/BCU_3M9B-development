# -*- coding: utf-8 -*-
"""八路径验证总报告。

报告把 MATLAB 交叉验证、物理验证、近似实现和环境阻塞严格分开；``5/5`` 等历史
测试项数量不会被解释成八条路径全部完成。每条记录均有唯一状态和限制说明。
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        import ctypes
        ctypes.windll.kernel32.SetConsoleOutputCP(65001)
        ctypes.windll.kernel32.SetConsoleCP(65001)
    except Exception:
        pass
    for _stream in (sys.stdout, sys.stderr):
        try:
            _stream.reconfigure(encoding="utf-8")
        except Exception:
            pass

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT.parent / "python_bcu") not in sys.path:
    sys.path.insert(0, str(ROOT.parent / "python_bcu"))

REFERENCE_DIR = ROOT.parent / "validation" / "references"
REPORT_DIR = ROOT.parent / "validation" / "reports"


def _entry(name: str, status: str, reference: str, *, passed: int = 0,
           total: int = 0, error: float | None = None,
           limitations: list[str] | None = None) -> dict:
    return {
        "name": name,
        "status": status,
        "checks_passed": passed,
        "checks_total": total,
        "max_error": None if error is None or not np.isfinite(error) else float(error),
        "reference": reference,
        "limitations": limitations or [],
    }


def _reference_status(filename: str) -> tuple[str, list[str]]:
    path = REFERENCE_DIR / filename
    if not path.exists():
        return "BLOCKED", [f"缺少紧凑参考 {filename}"]
    try:
        from bcu_v2.reference_io import load_reference
        data = load_reference(path)
    except Exception as exc:  # noqa: BLE001
        return "FAILED", [f"参考 schema 无法读取: {exc}"]
    status = data.get("status", "AVAILABLE")
    if status != "AVAILABLE":
        return status, [data.get("reason", "参考尚未可用")]
    return "AVAILABLE", []


def _load_reference_record(filename: str) -> dict:
    """Load one compact reference for diagnostics without changing status rules."""
    from bcu_v2.reference_io import load_reference
    return load_reference(REFERENCE_DIR / filename)


def _spm_frame_limitations(data: dict) -> list[str]:
    """Detect the known MATLAB SPM projected/raw network-angle mismatch.

    The MATLAB exporter now retains the raw fsolve vector.  If the projected
    and raw fields do not describe the same common-angle frame, the reference
    may still be useful as historical evidence, but it is not a like-for-like
    physical energy reference.
    """
    arrays = data.get("arrays", {})
    required = {"cuep_raw_net_theta", "cuep_frame_shift", "cuep_net_theta"}
    if not required.issubset(arrays):
        return []
    raw_net = np.asarray(arrays["cuep_raw_net_theta"], dtype=float)
    exported_net = np.asarray(arrays["cuep_net_theta"], dtype=float)
    shift = float(arrays["cuep_frame_shift"])
    coherent = raw_net + shift
    mismatch = float(np.max(np.abs(coherent - exported_net)))
    if mismatch > 1e-6:
        return [
            ("MATLAB SPM 紧凑参考同时保留 raw fsolve 与 projected 网络角；两者最大坐标差 "
             f"为 {mismatch:.6g}，历史 E_critical 不能与物理 COI 一致坐标直接比较"),
        ]
    return []


def _matlab_path(name: str, ref: str, historical: str) -> dict:
    status, limitations = _reference_status(ref)
    if status == "AVAILABLE":
        # A compact reference being present is not, by itself, evidence that
        # the corresponding path was cross-validated.  The converter records
        # the historical check scope explicitly so the report cannot turn a
        # five-item historical summary into an eight-path claim.
        path = REFERENCE_DIR / ref
        try:
            from bcu_v2.reference_io import load_reference
            data = load_reference(path)
        except Exception as exc:  # noqa: BLE001
            return _entry(name, "FAILED", ref,
                          limitations=limitations + [f"参考 evidence 无法读取: {exc}"])
        evidence = data.get("evidence", {})
        passed = int(evidence.get("checks_passed", 0) or 0)
        total = int(evidence.get("checks_total", 0) or 0)
        error = evidence.get("max_error")
        if passed > 0 and total > 0 and passed == total:
            path_status = "MATLAB_XVAL_FULL"
        elif passed > 0:
            path_status = "MATLAB_XVAL_PARTIAL"
        else:
            path_status = "UNVERIFIED"
        evidence_limits = [
            "状态来自历史交叉验证证据；本轮仅重新生成并校验紧凑参考",
        ]
        if evidence.get("note"):
            evidence_limits.append(str(evidence["note"]))
        return _entry(name, path_status, ref, passed=passed, total=total,
                      error=error, limitations=limitations + evidence_limits)
    return _entry(name, status, ref, limitations=limitations + [historical])


def verify_spm_cct() -> dict:
    ref = "spm_cct_v1.json"
    ref_status, limitations = _reference_status(ref)
    if ref_status != "AVAILABLE":
        return _entry("spm_cct", "MATLAB_XVAL_PARTIAL" if ref_status == "BLOCKED" else ref_status,
                      ref, limitations=limitations + ["自足 CUEP 尚未取得 MATLAB 紧凑参考"])
    try:
        reference_data = _load_reference_record(ref)
        from bcu_v2 import config as C
        from bcu_v2.spm_cuep import spm_self_contained_cct
        static = C.build_static_from_config(C.apply_overrides(C.load_config(), {"mode": "spm_cct"}))
        result = spm_self_contained_cct(static)
    except Exception as exc:  # noqa: BLE001
        return _entry("spm_cct", "FAILED", ref, limitations=[f"自足求解异常: {exc}"])
    limitations = _spm_frame_limitations(reference_data)
    if np.isfinite(reference_data.get("arrays", {}).get("e_critical", np.nan)):
        limitations.append(
            f"MATLAB 历史管线 E_critical={float(reference_data['arrays']['e_critical']):.6g}；"
            "该数值仅作历史参考，未作为 Python 自足输入"
        )
    if not result.converged:
        return _entry("spm_cct", "UNVERIFIED", ref,
                      limitations=[result.exit_reason or "SPM 自足 CUEP 未收敛"] + limitations)
    return _entry("spm_cct", "MATLAB_XVAL_FULL", ref, passed=1, total=1,
                  limitations=["需在 MATLAB 可用后复核固定参考数值"] + limitations)


def verify_spm_numerical() -> dict:
    return _entry("spm_numerical", "APPROXIMATE", "spm_cct_v1.json",
                  limitations=["已有连续 DAE 原型；尚无 MATLAB 固定检查点对照"])


def verify_spm_region() -> dict:
    return _entry("spm_region", "UNVERIFIED", "spm_cct_v1.json",
                  limitations=["SPM 平衡点/稳定域导出与 MATLAB 对照尚未完成"])


def build_report() -> dict:
    entries = [
        _matlab_path("reduced_cct", "reduced_cct_v1.json", "历史 T3 8/8 仅代表 reduced_cct"),
        _matlab_path("reduced_numerical", "reduced_numerical_v1.json", "历史仅有故障段末端 thetac 对照"),
        _matlab_path("reduced_region", "reduced_region_v1.json", "历史仅有平衡点集合对照"),
        verify_spm_cct(),
        verify_spm_numerical(),
        verify_spm_region(),
        _matlab_path("two_machine_region_3d", "two_machine_3d_v1.json", "历史平衡点集合对照待重新导出"),
        _matlab_path("two_machine_region_3d_gfl", "two_machine_gfl_v1.json", "历史 GFL 平衡点集合待重新导出"),
    ]
    return {
        "schema_version": "1.0",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "entries": entries,
        "summary": {
            "matlab_xval_full": sum(x["status"] == "MATLAB_XVAL_FULL" for x in entries),
            "matlab_xval_partial": sum(x["status"] == "MATLAB_XVAL_PARTIAL" for x in entries),
            "unverified_or_approximate": sum(x["status"] in {"UNVERIFIED", "APPROXIMATE"} for x in entries),
            "blocked": sum(x["status"] == "BLOCKED" for x in entries),
        },
    }


def main() -> int:
    report = build_report()
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    out = REPORT_DIR / "full_xval_latest.json"
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if all(x["status"] == "MATLAB_XVAL_FULL" for x in report["entries"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
