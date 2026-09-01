"""SPM 自足 CUEP 的先行测试。

测试只允许在 v2 内部生成临界能量；如果求解失败必须返回结构化失败状态，不能
从 MATLAB 参考或硬编码数字回填。
"""

from __future__ import annotations

import unittest
from unittest.mock import patch

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
from bcu_v2.spm_energy import (
    solve_spm_network,
    solve_spm_network_newton,
    spm_potential_energy,
)
from bcu_v2.spm_cuep import estimate_spm_fault_energy_peak


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
        if np.isfinite(result.omega_coi) and result.delta_gen.size:
            theta = result.theta_net
            voltage = result.voltage_net
            yfull = np.asarray(self.static.postfault.metadata["yfull_mod"])
            from bcu_v2.spm_energy import spm_generator_power
            pe = spm_generator_power(result.delta_gen, theta, voltage,
                                     yfull, self.static.preset.epu)
            expected = float(np.sum(self.static.preset.pmpu - pe) /
                             np.sum(self.static.preset.d))
            self.assertAlmostEqual(result.omega_coi, expected, places=8)
        if result.converged:
            self.assertTrue(np.isfinite(result.e_critical))
            self.assertLess(result.network_residual, 1e-8)

    def test_self_contained_result_supports_tuple_unpacking(self):
        result = spm_self_contained_cct(
            self.static, tfault=0.05, tunit=0.005, max_segments=1,
        )
        self.assertIsInstance(result, SpmSelfContainedResult)
        cuep, cct, ok = result
        self.assertIs(cuep, result.cuep)
        if np.isnan(result.cct):
            self.assertTrue(np.isnan(cct))
        else:
            self.assertEqual(cct, result.cct)
        self.assertEqual(ok, result.converged)
        self.assertFalse(result.used_external_ecritical)

    def test_zero_voltage_algebraic_root_is_not_physical_success(self):
        yfull = np.asarray(self.static.postfault.metadata["yfull_mod"])
        zero_guess = np.zeros(2 * (yfull.shape[0] - self.static.preset.ngen))
        _, converged, residual = solve_spm_network(
            self.static.postfault.sep_delta, yfull, self.static.preset.epu,
            guess=zero_guess,
        )
        self.assertFalse(converged)
        self.assertLess(residual, 1e-8)

    def test_matlab_newton_network_corrector_returns_physical_sep(self):
        yfull = np.asarray(self.static.postfault.metadata["yfull_mod"])
        z, converged, residual = solve_spm_network_newton(
            self.static.postfault.sep_delta,
            yfull,
            self.static.preset.epu,
            tol=1e-12,
        )
        self.assertTrue(converged, msg=f"Newton residual={residual:g}")
        self.assertLess(residual, 1e-10)
        self.assertTrue(np.all(z[6:] > 1e-4))

    def test_fault_energy_peak_uses_spm_dae_trajectory(self):
        """SPM energy must not silently fall back to the reduced trajectory."""
        with patch("bcu_3m9b.dynamics.integrate_reduced",
                   side_effect=AssertionError("reduced trajectory is not SPM evidence")):
            peak = estimate_spm_fault_energy_peak(
                self.static, tfault=0.05, tunit=0.005, max_points=16,
            )
        self.assertTrue(np.isfinite(peak))

    def test_fault_energy_reconstructs_postfault_network_states(self):
        """Energy peak must use MATLAB's postfault algebraic reconstruction."""
        peak = estimate_spm_fault_energy_peak(
            self.static, tfault=0.5, tunit=1e-3, max_points=64,
        )
        self.assertGreater(peak, 5.0)
        self.assertAlmostEqual(peak, 5.5678, delta=0.1)

    def test_path_energy_cal_reproduces_nonzero_ep4_without_changing_cuep_default(self):
        """MGP ray integration may request MATLAB's Ep4; CUEP default stays zero."""
        yfull = np.asarray(self.static.postfault.metadata["yfull_mod"])
        sep = np.asarray(self.static.postfault.sep_delta, dtype=float)
        sep_z, ok, residual = solve_spm_network(
            sep, yfull, self.static.preset.epu, tol=1e-12,
        )
        self.assertTrue(ok, msg=f"SEP residual={residual:g}")
        end = sep + np.array([0.05, -0.03, 0.08])
        end_z, ok, residual = solve_spm_network(
            end, yfull, self.static.preset.epu, guess=sep_z, tol=1e-12,
        )
        self.assertTrue(ok, msg=f"end residual={residual:g}")
        ep_default = spm_potential_energy(
            self.static.preset, self.static.postfault, yfull,
            sep, sep_z[:6], sep_z[6:], end, end_z[:6], end_z[6:],
        )
        ep_ray = spm_potential_energy(
            self.static.preset, self.static.postfault, yfull,
            sep, sep_z[:6], sep_z[6:], end, end_z[:6], end_z[6:],
            path_energy_cal=20,
        )
        self.assertAlmostEqual(ep_default[3], 0.0, places=12)
        self.assertTrue(np.isfinite(ep_ray[3]))
        self.assertGreater(abs(float(ep_ray[3])), 1e-8)


if __name__ == "__main__":
    unittest.main()
