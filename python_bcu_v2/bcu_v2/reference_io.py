"""读取和校验紧凑 MATLAB/Python 参考数据。

新格式是 JSON：顶层包含 ``schema_version``、``metadata`` 和 ``arrays``；复数数组使用
``{"real": ..., "imag": ...}``，不允许错误数据静默降级为零值。旧的 MATLAB ``.mat``
文件仍可读取，供过渡期诊断使用，但不会被新验证报告当作紧凑参考版本。
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np


def _is_complex_payload(value: Any) -> bool:
    return isinstance(value, dict) and "real" in value and "imag" in value


def _check_array(name: str, value: Any, errors: list[str]) -> None:
    if _is_complex_payload(value):
        real = np.asarray(value["real"])
        imag = np.asarray(value["imag"])
        if real.shape != imag.shape:
            errors.append(f"{name}: real/imag shape mismatch {real.shape} != {imag.shape}")
        return
    if isinstance(value, (list, tuple, int, float, bool)):
        try:
            np.asarray(value)
        except Exception as exc:  # noqa: BLE001
            errors.append(f"{name}: invalid array ({exc})")
        return
    errors.append(f"{name}: expected JSON array/scalar or complex payload")


def validate_reference_schema(data: dict[str, Any]) -> tuple[bool, list[str]]:
    """校验紧凑参考 schema，返回 ``(valid, errors)``。"""

    errors: list[str] = []
    if not isinstance(data, dict):
        return False, ["top-level value must be an object"]
    if data.get("schema_version") != "1.0":
        errors.append("schema_version must be '1.0'")
    metadata = data.get("metadata")
    if not isinstance(metadata, dict):
        errors.append("metadata must be an object")
    else:
        for key in ("case", "fault", "source_matlab_commit", "created_at"):
            if key not in metadata:
                errors.append(f"metadata.{key} is required")
    arrays = data.get("arrays")
    if not isinstance(arrays, dict):
        errors.append("arrays must be an object")
    else:
        for name, value in arrays.items():
            _check_array(f"arrays.{name}", value, errors)
    status = data.get("status", "AVAILABLE")
    if status not in {"AVAILABLE", "BLOCKED", "UNVERIFIED", "NOT_COMPARABLE"}:
        errors.append(f"unsupported status: {status}")
    if status != "AVAILABLE" and not data.get("reason"):
        errors.append("non-available reference must include reason")
    return not errors, errors


def load_reference(path: Path) -> dict[str, Any]:
    """加载 JSON 或旧 MATLAB ``.mat`` 参考；错误数据直接抛出异常。"""

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".json":
        data = json.loads(path.read_text(encoding="utf-8"))
        valid, errors = validate_reference_schema(data)
        if not valid:
            raise ValueError(f"invalid reference schema in {path}: {'; '.join(errors)}")
        return data
    if path.suffix.lower() == ".mat":
        try:
            from scipy.io import loadmat
        except ImportError as exc:  # pragma: no cover - environment dependent
            raise RuntimeError("scipy is required to read legacy .mat references") from exc
        raw = loadmat(str(path), squeeze_me=True, struct_as_record=False)
        return {key: value for key, value in raw.items() if not key.startswith("__")}
    raise ValueError(f"unsupported reference format: {path.suffix}")


def as_numpy(value: Any) -> np.ndarray:
    """把 schema 数组转为 ndarray；复数 payload 恢复为 complex 数组。"""

    if _is_complex_payload(value):
        return np.asarray(value["real"], dtype=float) + 1j * np.asarray(value["imag"], dtype=float)
    return np.asarray(value)
