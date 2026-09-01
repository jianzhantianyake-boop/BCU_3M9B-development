"""Convert legacy MATLAB verification snapshots into compact JSON references.

The source directory is read-only from the repository's point of view.  This
utility copies only selected arrays from legacy ``.mat`` files and records the
source file hash, source MATLAB commit and conversion metadata.  It never
copies a complete MATLAB workspace or a large trajectory into the repository.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat


SOURCE_COMMIT = "035f1475fd92e5639ff9b7fb78eb678ed2976e1c"
SCHEMA_VERSION = "1.0"

HISTORICAL_EVIDENCE: dict[str, dict[str, Any]] = {
    "reduced_cct_v1": {
        "checks_passed": 8,
        "checks_total": 8,
        "max_error": 1.0e-10,
        "kind": "historical_cross_validation",
        "note": "历史 T3 逐层对照摘要；本次只重新生成紧凑参考",
    },
    "reduced_numerical_v1": {
        "checks_passed": 1,
        "checks_total": 1,
        "max_error": 1.2e-11,
        "kind": "historical_cross_validation",
        "note": "故障段末端 thetac 代表性检查点",
    },
    "reduced_region_v1": {
        "checks_passed": 1,
        "checks_total": 1,
        "max_error": 4.5e-11,
        "kind": "historical_cross_validation",
        "note": "平衡点集合代表性对照；稳定边界曲线未覆盖",
    },
    "two_machine_3d_v1": {
        "checks_passed": 1,
        "checks_total": 1,
        "max_error": 6.2e-13,
        "kind": "historical_cross_validation",
        "note": "平衡点集合对照，D2 参数按 0.5",
    },
    "two_machine_gfl_v1": {
        "checks_passed": 1,
        "checks_total": 1,
        "max_error": 1.4e-11,
        "kind": "historical_cross_validation",
        "note": "GFL 平衡点集合对照",
    },
}


def _mat_struct(value: Any, name: str) -> Any:
    if not hasattr(value, name):
        raise KeyError(f"MATLAB field missing: {name}")
    return getattr(value, name)


def _array(value: Any) -> list[Any]:
    """Return JSON-compatible numeric arrays, preserving complex values."""
    a = np.asarray(value)
    if np.iscomplexobj(a):
        return {"real": np.real(a).tolist(), "imag": np.imag(a).tolist()}
    return a.tolist()


def _scalar(value: Any) -> Any:
    a = np.asarray(value)
    if a.size != 1:
        return _array(value)
    item = a.reshape(-1)[0]
    if isinstance(item, np.generic):
        item = item.item()
    return item


def _load(path: Path) -> dict[str, Any]:
    return loadmat(path, squeeze_me=True, struct_as_record=False)


def _reduced_cct(path: Path) -> dict[str, Any]:
    raw = _load(path)
    pre = raw["prefault"]
    post = raw["postfault"]
    fault = raw["fault"]
    critical = raw["Critical"]
    lea = _mat_struct(critical, "LEA")
    rea = _mat_struct(critical, "REA")
    return {
        "prefault_yred": _array(_mat_struct(pre, "Yred")),
        "fault_yred": _array(_mat_struct(fault, "Yred")),
        "postfault_yred": _array(_mat_struct(post, "Yred")),
        "prefault_sep_delta": _array(_mat_struct(pre, "SEP_delta")),
        "postfault_sep_delta": _array(_mat_struct(post, "SEP_delta")),
        "postfault_cuep_delta": _array(_mat_struct(post, "CUEP_delta")),
        "postfault_sep_residual": _array(_mat_struct(post, "SEP_Perr")),
        "postfault_cuep_residual": _array(_mat_struct(post, "CUEP_Perr")),
        "lea_cct": _scalar(_mat_struct(lea, "CCT")),
        "rea_cct": _scalar(_mat_struct(rea, "CCT")),
        "lea_exit_thetac": _array(_mat_struct(lea, "Exit_thetac")),
        "rea_exit_thetac": _array(_mat_struct(rea, "Exit_thetac")),
        "faultline": _array(_mat_struct(fault, "faultline")),
    }


def _reduced_numerical(path: Path) -> dict[str, Any]:
    raw = _load(path)
    return {name: _array(raw[name]) for name in ("theta_end", "omega_end", "thetac_end")}


def _reduced_region(path: Path) -> dict[str, Any]:
    raw = _load(path)
    return {"xeps": _array(raw["xeps"]), "flags": _array(raw["flags"])}


def _two_machine(path: Path) -> dict[str, Any]:
    raw = _load(path)
    xeps = raw["xeps"]
    # MATLAB cell arrays become object arrays under scipy; normalize each row.
    rows = [_array(row) for row in np.atleast_1d(xeps)]
    return {"xeps": rows, "flags": _array(raw["flags"])}


def _git_head(repo_root: Path) -> str:
    try:
        result = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "working-tree"


def _write_reference(
    output: Path,
    name: str,
    source: Path,
    arrays: dict[str, Any],
    *,
    repo_root: Path,
    created_at: str,
) -> dict[str, Any]:
    record = {
        "name": name,
        "schema_version": SCHEMA_VERSION,
        "status": "AVAILABLE",
        "reason": "从只读 MATLAB 紧凑快照提取；未保存完整工作区",
        "metadata": {
            "case": "case9_v2" if "two_machine" not in name else "two_machine_reference",
            "fault": "F9" if "two_machine" not in name else "not_applicable",
            "matlab_version": "R2024a (source snapshot metadata)",
            "source_matlab_commit": SOURCE_COMMIT,
            "source_file": source.name,
            "source_file_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "exporter_commit": _git_head(repo_root),
            "created_at": created_at,
        },
        "evidence": HISTORICAL_EVIDENCE.get(name, {}),
        "arrays": arrays,
    }
    output.write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
    return record


def convert(source_verify_root: Path, output_root: Path, repo_root: Path, *, replace_unverified: bool = False) -> dict[str, Any]:
    source_verify_root = source_verify_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    created_at = datetime.now(timezone.utc).isoformat()
    specs = [
        ("baseline_reduced.mat", "reduced_cct_v1.json", _reduced_cct),
        ("baseline_numerical.mat", "reduced_numerical_v1.json", _reduced_numerical),
        ("baseline_region.mat", "reduced_region_v1.json", _reduced_region),
        ("baseline_twomachine.mat", "two_machine_3d_v1.json", _two_machine),
        ("baseline_twomachine_gfl.mat", "two_machine_gfl_v1.json", _two_machine),
    ]
    entries: list[dict[str, Any]] = []
    for source_name, output_name, loader in specs:
        source = source_verify_root / source_name
        if not source.is_file():
            raise FileNotFoundError(source)
        output = output_root / output_name
        # A conversion is explicit and reproducible, but never silently overwrites
        # a reference that may have been generated by a newer exporter.
        if output.exists():
            if not replace_unverified:
                raise FileExistsError(f"refusing to overwrite existing reference: {output}")
            try:
                existing = json.loads(output.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise FileExistsError(f"existing reference is not a replaceable placeholder: {output}") from exc
            existing_source_hash = (
                existing.get("metadata", {}).get("source_file_sha256")
            )
            current_source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
            same_source_refresh = (
                existing.get("status") == "AVAILABLE"
                and existing_source_hash == current_source_hash
            )
            if existing.get("status") != "UNVERIFIED" and not same_source_refresh:
                raise FileExistsError(f"refusing to overwrite non-placeholder reference: {output}")
        record = _write_reference(output, output_name.removesuffix(".json"), source, loader(source), repo_root=repo_root, created_at=created_at)
        entries.append(
            {
                "name": record["name"],
                "path": output.name,
                "status": record["status"],
                "reference_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
                "source_file": source.name,
                "source_file_sha256": record["metadata"]["source_file_sha256"],
            }
        )
    # Preserve a separately exported current SPM reference when present.
    spm_path = output_root / "spm_cct_v1.json"
    if spm_path.is_file():
        spm = json.loads(spm_path.read_text(encoding="utf-8"))
        entries.append(
            {
                "name": spm.get("name", "spm_cct_v1"),
                "path": spm_path.name,
                "status": spm.get("status", "UNVERIFIED"),
                "reference_sha256": hashlib.sha256(spm_path.read_bytes()).hexdigest(),
                "source_file": "MATLAB batch export",
                "source_file_sha256": None,
            }
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "case": "case9_v2",
        "fault": "F9",
        "matlab_version": "R2024a (source snapshot metadata)",
        "exporter_commit": _git_head(repo_root),
        "source_matlab_commit": SOURCE_COMMIT,
        "created_at": created_at,
        "entries": entries,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-verify-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--replace-unverified", action="store_true")
    args = parser.parse_args()
    manifest = convert(args.source_verify_root, args.output_root, args.repo_root, replace_unverified=args.replace_unverified)
    manifest_path = args.output_root / "reference_manifest.json"
    if manifest_path.exists() and not args.replace_unverified:
        raise FileExistsError(f"refusing to overwrite existing manifest: {manifest_path}")
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
