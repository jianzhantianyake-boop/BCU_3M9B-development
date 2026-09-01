"""关键方程变异测试：每个受控错误都必须被至少一个不变量检测。"""

from __future__ import annotations

import unittest

import numpy as np

import bcu_v2  # noqa: F401
from bcu_3m9b import build_static_result
from bcu_3m9b.dynamics import reduced_rhs
from bcu_3m9b.equilibrium import electrical_power
from bcu_3m9b.network import remove_fault_line
from bcu_v2 import spm_energy
from bcu_v2.spm_region import branch_continuity


class MutationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.static = build_static_result()
        cls.preset = cls.static.preset

    def test_kron_sign_mutation(self):
        state = self.static.postfault
        bad = state.ynn + state.ynr @ np.linalg.solve(state.yrr, state.yrn)
        self.assertGreater(np.max(np.abs(bad - state.yred)), 1e-6)

    def test_swing_acceleration_sign_mutation(self):
        s = self.static
        x = np.r_[s.prefault.sep_delta, np.full(3, s.basevalue.omega_b)]
        good = reduced_rhs(x, s.prefault.yred, self.preset, s.basevalue)
        bad = np.r_[good[:3], -good[3:]]
        self.assertGreater(np.max(np.abs(good - bad)), 1e-6)

    def test_fault_line_removal_mutation(self):
        case = s = self.static.case
        good = remove_fault_line(case.branch, self.preset.fault_line)
        self.assertLess(good.shape[0], case.branch.shape[0])

    def test_coi_weight_mutation(self):
        theta = np.array([0.2, 0.7, -0.1])
        good = theta - np.dot(self.preset.m, theta) / np.sum(self.preset.m)
        bad = theta - np.mean(theta)
        self.assertGreater(np.max(np.abs(good - bad)), 1e-4)

    def test_generator_order_mutation(self):
        s = self.static
        theta = s.postfault.sep_delta
        good = electrical_power(theta, s.postfault.yred, self.preset.epu)
        bad = electrical_power(theta, s.postfault.yred, self.preset.epu[[1, 0, 2]])
        self.assertGreater(np.max(np.abs(good - bad)), 1e-5)

    def test_energy_angle_sign_mutation(self):
        s = self.static
        sep = s.postfault.sep_delta
        end = sep + np.array([0.1, -0.2, 0.1])
        good = np.sum(__import__("bcu_3m9b.energy", fromlist=["potential_energy"]).potential_energy(
            self.preset, s.postfault, sep, end))
        self.assertGreater(abs(good - (-good)), 1e-8)

    def test_spm_pq_sign_mutation(self):
        s = self.static
        y = np.asarray(s.postfault.metadata["yfull_mod"])
        dg = s.postfault.sep_delta
        x, ok, _ = spm_energy.solve_spm_network(dg, y, self.preset.epu)
        self.assertTrue(ok)
        x_bad = x.copy(); x_bad[0] += 0.1
        residual = spm_energy.spm_network_residual(x_bad, dg, y, self.preset.epu)
        self.assertGreater(np.max(np.abs(residual - (-residual))), 1e-10)

    def test_voltage_angle_mutation(self):
        s = self.static
        y = np.asarray(s.postfault.metadata["yfull_mod"])
        dg = s.postfault.sep_delta
        x, ok, _ = spm_energy.solve_spm_network(dg, y, self.preset.epu)
        self.assertTrue(ok)
        wrong = np.r_[x[:6], np.ones(6)]
        self.assertGreater(np.linalg.norm(spm_energy.spm_network_residual(
            wrong, dg, y, self.preset.epu)), 1e-3)

    def test_branch_continuity_mutation(self):
        states = np.zeros((3, 4)); states[2, 0] = 10.0
        ok, jump = branch_continuity(states, tolerance=1.0)
        self.assertFalse(ok)
        self.assertGreater(jump, 1.0)


if __name__ == "__main__":
    unittest.main()
