"""reference_io 的 schema/错误处理测试。"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from bcu_v2.reference_io import as_numpy, load_reference, validate_reference_schema


def _valid_reference() -> dict:
    return {
        "schema_version": "1.0",
        "status": "AVAILABLE",
        "metadata": {
            "case": "case9_v2",
            "fault": "F9",
            "source_matlab_commit": "035f1475fd92e5639ff9b7fb78eb678ed2976e1c",
            "created_at": "2026-09-01T00:00:00Z",
        },
        "arrays": {
            "Yred": {"real": [[1.0, 0.0]], "imag": [[0.0, -1.0]]},
            "sep_delta": [0.1, 0.2, 0.3],
        },
    }


class ReferenceIoTests(unittest.TestCase):
    def test_valid_schema_and_complex_roundtrip(self) -> None:
        data = _valid_reference()
        self.assertEqual(validate_reference_schema(data), (True, []))
        value = as_numpy(data["arrays"]["Yred"])
        np.testing.assert_allclose(value, np.array([[1.0 + 0j, -1j]]))

    def test_invalid_schema_is_rejected(self) -> None:
        data = _valid_reference()
        data["arrays"]["Yred"]["imag"] = []
        valid, errors = validate_reference_schema(data)
        self.assertFalse(valid)
        self.assertTrue(any("shape mismatch" in error for error in errors))

    def test_load_json_does_not_silently_zero_bad_data(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.json"
            data = _valid_reference()
            data.pop("metadata")
            path.write_text(json.dumps(data), encoding="utf-8")
            with self.assertRaises(ValueError):
                load_reference(path)

    def test_blocked_reference_requires_reason(self) -> None:
        data = _valid_reference()
        data["status"] = "BLOCKED"
        valid, errors = validate_reference_schema(data)
        self.assertFalse(valid)
        self.assertIn("non-available reference must include reason", errors)

    def test_current_spm_reference_exposes_raw_frame_diagnostic(self) -> None:
        """The historical MATLAB energy reference must not hide a frame mismatch."""
        path = Path(__file__).resolve().parents[1] / "validation" / "references" / "spm_cct_v1.json"
        if not path.exists():
            self.skipTest("compact SPM reference has not been exported")
        data = load_reference(path)
        arrays = data["arrays"]
        required = {"cuep_raw_net_theta", "cuep_frame_shift", "cuep_net_theta"}
        self.assertTrue(required.issubset(arrays))
        raw = np.asarray(arrays["cuep_raw_net_theta"], dtype=float)
        projected = np.asarray(arrays["cuep_net_theta"], dtype=float)
        coherent = raw + float(arrays["cuep_frame_shift"])
        self.assertGreater(float(np.max(np.abs(coherent - projected))), 1e-6)


if __name__ == "__main__":
    unittest.main()
