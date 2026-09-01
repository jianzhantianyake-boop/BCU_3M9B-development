"""阻尼研究驱动的接口回归测试。"""

from __future__ import annotations

import sys
import unittest
from types import SimpleNamespace
from unittest.mock import patch

ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_bcu_v2"))
sys.path.insert(0, str(ROOT / "python_bcu"))

import run_damping_study as study
from bcu_v2 import config


class DampingStudyTests(unittest.TestCase):
    def test_batch_gate_blocks_without_creating_results(self):
        ok, reasons = study.gate_status()
        self.assertFalse(ok)
        self.assertTrue(reasons)
        before = study.OUT.exists()
        self.assertEqual(study.main(), 2)
        self.assertEqual(study.OUT.exists(), before)

    def test_case_runner_uses_self_contained_converged_flag(self):
        case = study.build_cases()[0]
        fake_static = object()
        fake_reduced = {"lea": SimpleNamespace(cct=0.2),
                        "rea_cct": 0.21}
        fake_uep = SimpleNamespace(found=True, v_cuep=1.2)
        fake_spm = SimpleNamespace(
            cuep=SimpleNamespace(e_critical=1.3, equilibrium_residual=1e-10,
                                 network_residual=2e-10, equilibrium_type="type-1"),
            cct=0.19,
            converged=True,
            exit_reason="",
        )
        with patch.object(config, "load_config", return_value={}), \
                patch.object(config, "apply_overrides", return_value={}), \
                patch.object(config, "build_static_from_config", return_value=fake_static), \
                patch.object(study, "_source_revision", return_value="test-revision"), \
                patch("bcu_v2.fixes.run_experiment_clean", return_value=fake_reduced), \
                patch("bcu_v2.cuep.controlling_uep", return_value=fake_uep), \
                patch("bcu_v2.spm_cuep.spm_self_contained_cct", return_value=fake_spm):
            row = study._run_case(case)
        self.assertEqual(row["convergence_status"], "PASSED")
        self.assertEqual(row["spm_lea_cct"], 0.19)
        self.assertIsNone(row["failure_reason"])


if __name__ == "__main__":
    unittest.main()
