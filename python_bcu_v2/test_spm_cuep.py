"""SPM 自足 CUEP 的先行测试。

测试只允许在 v2 内部生成临界能量；如果求解失败必须返回结构化失败状态，不能
从 MATLAB 参考或硬编码数字回填。
"""

from __future__ import annotations

import unittest

import numpy as np

import bcu_v2  # noqa: F401  # 触发 v1 sibling path 引导
from bcu_3m9b import build_static_result
from bcu_v2.spm_cuep import (
    SpmCuepResult,
    SpmMgpResult,
    SpmSelfContainedResult,
    solve_spm_cuep,
    spm_self_contained_cct,
    trace_spm_mgp,
)


class SpmCuepTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.static = build_static_result()

    def test_trace_returns_structured_state(self):
        result = trace_spm_mgp(self.static, max_segments=2, segment_steps=2)
        self.assertIsInstance(result, SpmMgpResult)
        self.assertIsInstance(result.converged, bool)
        self.assertTrue(result.exit_reason)

    def test_solver_has_no_external_ecritical_dependency(self):
        mgp = trace_spm_mgp(self.static, max_segments=4, segment_steps=2)
        result = solve_spm_cuep(self.static, mgp)
        self.assertIsInstance(result, SpmCuepResult)
        self.assertFalse(result.used_external_ecritical)
        if result.converged:
            self.assertTrue(np.isfinite(result.e_critical))
            self.assertLess(result.network_residual, 1e-8)

    def test_self_contained_result_supports_tuple_unpacking(self):
        result = spm_self_contained_cct(self.static, tfault=0.05, tunit=0.005)
        self.assertIsInstance(result, SpmSelfContainedResult)
        cuep, cct, ok = result
        self.assertIs(cuep, result.cuep)
        if np.isnan(result.cct):
            self.assertTrue(np.isnan(cct))
        else:
            self.assertEqual(cct, result.cct)
        self.assertEqual(ok, result.converged)
        self.assertFalse(result.used_external_ecritical)


if __name__ == "__main__":
    unittest.main()
