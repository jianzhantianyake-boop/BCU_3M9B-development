"""SPM 平衡点输出结构测试。"""

from __future__ import annotations

import unittest

import numpy as np
import bcu_v2  # noqa: F401
from bcu_3m9b import build_static_result
from bcu_v2.spm_region import (SpmEquilibrium, branch_continuity,
                               enumerate_spm_equilibria,
                               trace_spm_stable_manifold)


class SpmRegionTests(unittest.TestCase):
    def test_equilibrium_records_include_branch_and_residual(self):
        static = build_static_result()
        records = enumerate_spm_equilibria(static, grid_points=5)
        self.assertTrue(records)
        self.assertTrue(all(isinstance(item, SpmEquilibrium) for item in records))
        self.assertTrue(all(item.branch_id and item.residual_norm < 1e-6 for item in records))

    def test_continuous_network_branch_is_registered(self):
        records = enumerate_spm_equilibria(build_static_result(), grid_points=21)
        self.assertTrue(records)
        type1 = [item for item in records if item.equilibrium_type == "type-1"]
        self.assertTrue(type1)
        self.assertLess(type1[0].continuity_error, 1.0)
        states = __import__("numpy").array([[0.0, 0.0], [0.1, -0.1], [0.2, -0.2]])
        ok, maximum = branch_continuity(states, tolerance=1.0)
        self.assertTrue(ok)
        self.assertGreater(maximum, 0.0)

    def test_stable_manifold_trace_keeps_algebraic_residual_small(self):
        static = build_static_result()
        type1 = [item for item in enumerate_spm_equilibria(static, grid_points=21)
                 if item.equilibrium_type == "type-1"]
        self.assertEqual(len(type1), 1)
        result = trace_spm_stable_manifold(
            static, type1[0], np.array([-0.4288739686, 0.9033643335]),
            -1, sample_times=np.array([0.0, 0.01]),
        )
        self.assertTrue(result["converged"])
        self.assertLess(float(np.max(result["residual_norm"])), 1e-6)


if __name__ == "__main__":
    unittest.main()
