"""SPM 完整联合根审计接口的先行测试。"""

from __future__ import annotations

import sys
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = __import__("pathlib").Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT.parent / "python_bcu") not in sys.path:
    sys.path.insert(0, str(ROOT.parent / "python_bcu"))

from bcu_3m9b import build_static_result
from run_spm_joint_audit import load_seed_file
from bcu_v2.spm_joint_audit import audit_spm_joint_roots


class SpmJointAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.static = build_static_result()

    def test_audit_reports_structured_counts_without_external_energy(self):
        report = audit_spm_joint_roots(self.static, max_starts=2, random_seed=20260902)
        self.assertEqual(report["starts_requested"], 2)
        self.assertEqual(report["starts_evaluated"], 2)
        self.assertFalse(report["used_external_ecritical"])
        self.assertIn("roots", report)
        for root in report["roots"]:
            self.assertTrue(np.isfinite(root["network_residual"]))
            self.assertTrue(np.isfinite(root["equilibrium_residual"]))
            self.assertTrue(np.isfinite(root["e_critical"]))

    def test_audit_preserves_continuous_angle_branch(self):
        report = audit_spm_joint_roots(self.static, max_starts=1, random_seed=20260902)
        self.assertIn("angle_wraps", report)
        self.assertIsInstance(report["angle_wraps"], list)

    def test_audit_accepts_explicit_physical_network_seed(self):
        # 该状态来自只读 MATLAB 诊断导出；它用于验证审计器可以重放一个
        # 已知的物理网络分支，而不是把随机初值搜索误报成穷尽搜索。
        state = {
            "delta_gen": [-0.7585584172568889, 1.8575840695226964, 1.9978354475151372],
            "theta_net": [
                -0.7126093610476225, -0.3030526483798677,
                -0.7780140476285538, 1.5615857062985627,
                1.6256814059614753, 1.7637103117045405,
            ],
            "voltage_net": [
                0.6115351306051858, 0.2912895920212790,
                0.5929088068644735, 0.5805594920242412,
                0.6514770252835123, 0.7739313040016673,
            ],
        }
        report = audit_spm_joint_roots(
            self.static, max_starts=1, seed_states=[state], random_seed=20260902,
        )
        self.assertEqual(report["starts_evaluated"], 1)
        self.assertEqual(report["unique_converged_root_count"], 1)
        root = report["roots"][0]
        self.assertEqual(root["equilibrium_type"], "type-1")
        self.assertLess(root["network_residual"], 1e-8)
        self.assertAlmostEqual(root["e_critical"], 7.5766832677, places=8)

    def test_seed_file_extracts_mgp_deltas_and_physical_state(self):
        payload = {
            "trace": [
                {"update_delta": [1, 2, 3], "end_delta": [4, 5, 6]},
                {"update_delta": [7, 8, 9]},
            ],
            "physical_cuep": {
                "cuep_delta": [10, 11, 12],
                "corrected_net_theta": [0, 0, 0, 0, 0, 0],
                "cuep_net_voltage": [1, 1, 1, 1, 1, 1],
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "seeds.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            deltas, states = load_seed_file(path)
        self.assertEqual(len(deltas), 3)
        self.assertEqual(len(states), 1)
        self.assertEqual(states[0]["delta_gen"], [10, 11, 12])


if __name__ == "__main__":
    unittest.main()
