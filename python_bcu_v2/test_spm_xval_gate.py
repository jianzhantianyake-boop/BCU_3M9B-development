"""SPM 固定检查点参考的可比性门禁测试。"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

if sys.platform == "win32":
    try:
        import ctypes
        ctypes.windll.kernel32.SetConsoleOutputCP(65001)
    except Exception:
        pass

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT.parent / "python_bcu") not in sys.path:
    sys.path.insert(0, str(ROOT.parent / "python_bcu"))

import run_full_xval


class SpmXvalGateTests(unittest.TestCase):
    def test_spm_cct_does_not_promote_historical_mixed_frame_energy(self):
        entry = run_full_xval.verify_spm_cct()
        self.assertEqual(entry["status"], "MATLAB_XVAL_PARTIAL")
        self.assertTrue(any("待导入同坐标 MATLAB 紧凑参考" in item
                            for item in entry["limitations"]))
        self.assertTrue(any("本轮 Python 自足候选 E_critical" in item
                            for item in entry["limitations"]))
        self.assertTrue(all("3.375" not in item or "历史" in item
                            for item in entry["limitations"]))

    def test_matlab_cuep_reference_must_pass_physical_network_gate(self):
        diagnostics = run_full_xval.inspect_spm_cuep_reference(
            ROOT.parent / "validation" / "references" / "spm_cct_v1.json"
        )
        self.assertFalse(diagnostics["comparable"])
        self.assertGreater(diagnostics["network_residual"], 1e-6)
        entry = run_full_xval.verify_spm_cct()
        self.assertTrue(any("network residual" in item for item in entry["limitations"]))

    def test_matlab_cuep_raw_frame_correction_is_physical(self):
        diagnostics = run_full_xval.inspect_spm_cuep_reference(
            ROOT.parent / "validation" / "references" / "spm_cct_v1.json"
        )
        self.assertIn("corrected_raw_network_residual", diagnostics)
        self.assertLess(diagnostics["corrected_raw_network_residual"], 1e-6)
        self.assertIn("corrected_raw_e_critical", diagnostics)
        self.assertAlmostEqual(diagnostics["corrected_raw_e_critical"], 7.57668, delta=1e-2)

    def test_spm_cct_energy_gate_uses_matlab_half_second_window(self):
        entry = run_full_xval.verify_spm_cct()
        self.assertTrue(any("本轮故障能量峰值=5.567" in item
                            for item in entry["limitations"]))

    def test_matlab_fault_reference_is_flagged_when_network_residual_is_large(self):
        path = ROOT.parent / "validation" / "references" / "spm_numerical_v1.json"
        diagnostics = run_full_xval.inspect_spm_fault_reference(path)
        self.assertFalse(diagnostics["comparable"])
        self.assertGreater(diagnostics["max_fault_residual"], 1e-6)
        self.assertIn("fault1", diagnostics["reason"])

    def test_fault1_v2_reference_is_comparable_and_python_matches(self):
        path = ROOT.parent / "validation" / "references" / "spm_numerical_v2.json"
        if not path.exists():
            self.skipTest("fault1 compact SPM reference has not been exported")
        diagnostics = run_full_xval.inspect_spm_fault_reference(path)
        self.assertTrue(diagnostics["comparable"])
        self.assertLess(diagnostics["max_fault_residual"], 1e-6)
        entry = run_full_xval.verify_spm_numerical()
        self.assertEqual(entry["status"], "MATLAB_XVAL_FULL")
        self.assertEqual(entry["checks_passed"], entry["checks_total"])
        self.assertLess(entry["max_error"], 1e-6)

    def test_spm_region_compares_equilibria_and_manifold_checkpoints(self):
        entry = run_full_xval.verify_spm_region()
        self.assertEqual(entry["status"], "MATLAB_XVAL_FULL")
        self.assertGreaterEqual(entry["checks_total"], 20)
        self.assertEqual(entry["checks_passed"], entry["checks_total"])
        self.assertLess(entry["max_error"], 1e-6)
        self.assertTrue(any("固定采样点" in item for item in entry["limitations"]))


if __name__ == "__main__":
    unittest.main()
