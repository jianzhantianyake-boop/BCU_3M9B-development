"""SPM 联合同伦诊断的先行测试。"""

from __future__ import annotations

import sys
import unittest

import numpy as np

ROOT = __import__("pathlib").Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT.parent / "python_bcu") not in sys.path:
    sys.path.insert(0, str(ROOT.parent / "python_bcu"))

from bcu_3m9b import build_static_result
from bcu_v2.spm_homotopy import continue_spm_joint_homotopy


class SpmHomotopyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.static = build_static_result()

    def test_rejects_incomplete_or_nonfinite_start_state(self):
        with self.assertRaises(ValueError):
            continue_spm_joint_homotopy(self.static, {"delta_gen": [0, 1, 2]})
        with self.assertRaises(ValueError):
            continue_spm_joint_homotopy(
                self.static,
                {
                    "delta_gen": [0, 1, np.nan],
                    "theta_net": [0] * 6,
                    "voltage_net": [1] * 6,
                },
            )

    def test_mgp_seed_continues_to_physical_joint_root_without_external_energy(self):
        result = continue_spm_joint_homotopy(
            self.static,
            {
                "delta_gen": [-0.7512082360481194, 1.8382429209826538, 1.9813283429391368],
                "theta_net": [
                    -0.697282209694518, -0.2800030456728983,
                    -0.7626868962754495, 1.5335837610837106,
                    1.6002976966655214, 1.7416426781159957,
                ],
                "voltage_net": [
                    0.6158023969067533, 0.3015138536317254,
                    0.5970460994660146, 0.5848870287656607,
                    0.6544477867376334, 0.775593528922797,
                ],
            },
            steps=10,
        )
        self.assertTrue(result.converged)
        self.assertFalse(result.used_external_ecritical)
        self.assertLess(result.network_residual, 1e-8)
        self.assertLess(result.equilibrium_residual, 1e-8)
        self.assertTrue(np.all(result.voltage_net > 1e-4))
        self.assertAlmostEqual(result.e_critical, 7.5766832677, places=6)


if __name__ == "__main__":
    unittest.main()
