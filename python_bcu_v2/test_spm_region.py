"""SPM 平衡点输出结构测试。"""

from __future__ import annotations

import unittest

import bcu_v2  # noqa: F401
from bcu_3m9b import build_static_result
from bcu_v2.spm_region import SpmEquilibrium, enumerate_spm_equilibria


class SpmRegionTests(unittest.TestCase):
    def test_equilibrium_records_include_branch_and_residual(self):
        static = build_static_result()
        records = enumerate_spm_equilibria(static, grid_points=5)
        self.assertTrue(records)
        self.assertTrue(all(isinstance(item, SpmEquilibrium) for item in records))
        self.assertTrue(all(item.branch_id and item.residual_norm < 1e-6 for item in records))


if __name__ == "__main__":
    unittest.main()
