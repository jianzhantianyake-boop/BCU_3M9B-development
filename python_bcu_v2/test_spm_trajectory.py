"""SPM 四阶段轨迹与 RK45/Radau 一致性测试。"""

from __future__ import annotations

import unittest

import numpy as np

import bcu_v2  # noqa: F401
from bcu_3m9b import build_static_result
from bcu_v2.spm_dae import SpmTrajectoryResult, simulate_spm_trajectory


class SpmTrajectoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.static = build_static_result()

    def test_four_phase_result_and_residual(self):
        result = simulate_spm_trajectory(self.static, clear_time=0.03,
                                          postfault_time=0.03, tunit=0.005)
        self.assertIsInstance(result, SpmTrajectoryResult)
        self.assertEqual(result.time.shape[0], result.delta_gen.shape[0])
        self.assertEqual(result.time.shape[0], result.phase_labels.shape[0])
        self.assertEqual(result.phase_labels[0], "prefault")
        self.assertIn("clearing", set(result.phase_labels))
        finite = result.algebraic_residual[np.isfinite(result.algebraic_residual)]
        self.assertGreater(finite.size, 0)
        self.assertLess(float(np.max(finite)), 1e-7)

    def test_rk45_and_radau_are_close(self):
        rk = simulate_spm_trajectory(self.static, clear_time=0.02,
                                     postfault_time=0.02, tunit=0.005, method="RK45")
        rad = simulate_spm_trajectory(self.static, clear_time=0.02,
                                      postfault_time=0.02, tunit=0.005, method="Radau")
        self.assertLess(float(np.max(np.abs(rk.delta_gen - rad.delta_gen))), 2e-4)


if __name__ == "__main__":
    unittest.main()
