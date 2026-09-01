"""SPM 四阶段轨迹与 RK45/Radau 一致性测试。"""

from __future__ import annotations

import unittest

import numpy as np

import bcu_v2  # noqa: F401
from bcu_3m9b import build_static_result
from bcu_3m9b.spm import algebraic_residual
from bcu_v2.spm_dae import (SpmTrajectoryResult, select_spm_checkpoints,
                             _algebraic_context, _make_solver,
                             simulate_spm_trajectory)


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

    def test_fixed_checkpoints_do_not_reuse_out_of_range_state(self):
        clear_time = 0.03
        postfault_time = 0.03
        tunit = 0.005
        result = simulate_spm_trajectory(self.static, clear_time=clear_time,
                                         postfault_time=postfault_time, tunit=tunit)
        checkpoints = select_spm_checkpoints(
            result, clear_time=clear_time, postfault_time=postfault_time,
            tunit=tunit,
        )
        self.assertEqual([item["label"] for item in checkpoints], [
            "t0", "pre-clearing", "clearing", "post-clearing-10ms",
            "post-clearing-50ms", "post-clearing-100ms", "final",
        ])
        self.assertTrue(all(item["available"] for item in checkpoints[:4]))
        self.assertFalse(checkpoints[4]["available"])
        self.assertFalse(checkpoints[5]["available"])
        self.assertTrue(checkpoints[6]["available"])
        self.assertAlmostEqual(checkpoints[2]["actual_time"], clear_time, places=12)

    def test_fault_network_solver_rejects_inconsistent_cold_root(self):
        """SPM fault equations must use the loaded impedance network and zero PQ load."""
        preset = self.static.preset
        yfull, load_pq, ngen, nload = _algebraic_context(self.static.fault, preset)
        solver = _make_solver(yfull, load_pq, ngen, nload)
        z = solver(self.static.prefault.sep_delta)
        residual = np.linalg.norm(
            algebraic_residual(z, self.static.prefault.sep_delta,
                               yfull, load_pq, ngen)
        )
        self.assertLess(residual, 1e-7)
        self.assertTrue(np.all(np.isfinite(z[nload:])))
        self.assertTrue(np.all(z[nload:] > 1e-4))


if __name__ == "__main__":
    unittest.main()
