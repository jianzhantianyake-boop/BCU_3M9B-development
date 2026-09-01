# -*- coding: utf-8 -*-
"""P1.3: 严格 DAE 级的结构保持模型(SPM), 根治代数解脆弱.

思路(约束流形上的降阶 DAE):
    结构保持模型是 index-1 DAE: 发电机 (δg, ωg) 为微分变量, 负荷母线 (δL, V) 为代数变量,
    受功率平衡约束 g(δg, z)=0 约束. 这里把它当"约束流形上的降阶 ODE"求解:
      - 用 scipy.solve_ivp 自适应积分发电机状态(method='RK45' 或刚性 'Radau');
      - 每次 RHS 求值时用**连续法**(上一步解热启动的 scipy.root 校正)把 g=0 解到机器精度,
        使代数约束在每个求值点都严格满足. 相比 v1"每步一次冷启动牛顿", 连续法的初值始终贴近
        真解, 从而根治偶发不收敛.

对比 v1:
    v1.simulate_spm: 固定步 + 每步冷启动/上一步热启动牛顿, 偶发首步不收敛即抛异常.
    本模块: 自适应步 + 连续法热启动 + scipy.root + 失败回退(default guess / 同伦), 更稳更准.

依赖: scipy.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root


@dataclass
class SpmTrajectoryResult:
    """统一的 SPM 四阶段轨迹结果。"""

    time: np.ndarray
    delta_gen: np.ndarray
    omega_gen: np.ndarray
    theta_net: np.ndarray
    voltage_net: np.ndarray
    algebraic_residual: np.ndarray
    phase_labels: np.ndarray
    converged: bool


def _algebraic_context(state, preset):
    """内部: 取 MATLAB SPM 所用的重排恒阻抗网络上下文.

    MATLAB 的 SPM 在 ``Yfull_mod`` 中已经并入恒阻抗负荷，因此网络母线的
    P/Q 代数约束为零；它不是 v1 的 ``无负荷导纳 + 恒功率负荷`` 方程。故障
    工况还会删去故障母线，``yfull_mod`` 的宽度随工况变化，但发电机始终在
    前 ``ngen`` 个节点。
    """

    yfull_mod = state.metadata.get("yfull_mod")
    if yfull_mod is None:
        # Compatibility fallback for hand-built NetworkState values.  The
        # static builder always supplies yfull_mod, so this branch is only a
        # diagnostic convenience and does not silently add constant-power
        # loads.
        yfull_mod = state.metadata.get("yorg_mod", state.yfull)
    yfull = np.asarray(yfull_mod, dtype=complex)
    n = preset.ngen
    nload = yfull.shape[0] - n
    if nload < 0:
        raise ValueError("SPM network has fewer nodes than generators")
    load_pq = np.zeros((nload, 2), dtype=float)
    return yfull, load_pq, n, nload


def _make_solver(yorg, load_pq, ngen, nload, tol=1e-11, initial_guess=None,
                 epu: np.ndarray | None = None):
    """内部: 返回带连续法热启动 + 回退的代数求解闭包 solve(δg)->z."""

    from bcu_3m9b.spm import algebraic_residual
    from .spm_energy import spm_network_residual

    default_guess = np.r_[np.zeros(nload), np.ones(nload)]
    if initial_guess is not None:
        candidate = np.asarray(initial_guess, dtype=float).reshape(-1)
        if candidate.size != 2 * nload:
            raise ValueError("algebraic_guess has incompatible network width")
        default_guess = candidate.copy()
    cache = {"z": default_guess.copy()}

    def solve(delta_gen: np.ndarray) -> np.ndarray:
        dg = np.asarray(delta_gen, dtype=float)

        def resid(z):
            # The strict SPM DAE must use the generator internal-voltage
            # magnitudes Epu.  The legacy v1 helper assumes unit magnitudes;
            # retain it only for direct compatibility callers that omit epu.
            if epu is None:
                return algebraic_residual(z, dg, yorg, load_pq, ngen)
            return spm_network_residual(z, dg, yorg, np.asarray(epu, dtype=float), load_pq)

        # 校正: 先用上一步解(连续法热启动).
        sol = root(resid, cache["z"], method="hybr", tol=tol)
        z, r = sol.x, float(np.linalg.norm(resid(sol.x)))
        if r > 1e-7:  # 回退1: 默认初值
            s2 = root(resid, default_guess, method="hybr", tol=tol)
            if np.linalg.norm(resid(s2.x)) < r:
                z, r = s2.x, float(np.linalg.norm(resid(s2.x)))
        if r > 1e-7:  # 回退2: lm(阻尼最小二乘)
            s3 = root(resid, cache["z"], method="lm")
            if np.linalg.norm(resid(s3.x)) < r:
                z, r = s3.x, float(np.linalg.norm(resid(s3.x)))
        cache["z"] = z
        return z

    def reset(z=None):
        cache["z"] = default_guess.copy() if z is None else np.asarray(z, dtype=float)

    solve.reset = reset  # type: ignore[attr-defined]
    return solve


def _strict_spm_generator_rhs(delta_gen: np.ndarray, omega: np.ndarray,
                              state, preset, basevalue, algebraic: np.ndarray) -> np.ndarray:
    """MATLAB-compatible SPM generator RHS using explicit internal voltages."""

    yorg, _load_pq, ngen, nload = _algebraic_context(state, preset)
    delta_gen = np.asarray(delta_gen, dtype=float)
    omega = np.asarray(omega, dtype=float)
    z = np.asarray(algebraic, dtype=float)
    theta = np.r_[delta_gen, z[:nload]]
    voltage = np.r_[np.asarray(preset.epu, dtype=float), z[nload:]]
    phasor = voltage * np.exp(1j * theta)
    injection = phasor * np.conj(yorg @ phasor)
    pe = injection.real[:ngen]
    coi = float(np.dot(omega, preset.m) / np.sum(preset.m))
    return np.r_[omega - coi,
                 (np.asarray(preset.pmpu, dtype=float) - pe
                  - np.asarray(preset.d, dtype=float) * (omega - basevalue.omega_b))
                 / np.asarray(preset.m, dtype=float)]


def simulate_spm_dae(tlength: float, tunit: float, state, preset, basevalue,
                     delta0: np.ndarray, omega0: np.ndarray,
                     method: str = "RK45", rtol: float = 1e-8, atol: float = 1e-10,
                     algebraic_guess: np.ndarray | None = None) -> Dict:
    """严格 DAE 级 SPM 仿真(约束流形降阶 ODE + 连续法).

    使用方法:
        传入时长/输出步长/网络工况/参数/基值/初值; method 可选 'RK45'(默认, 快)或
        'Radau'/'BDF'(刚性). 返回 {time, delta, omega, delta_coi, algebraic, success, method}.
    """

    from bcu_3m9b.spm import spm_generator_rhs

    yorg, load_pq, n, nload = _algebraic_context(state, preset)
    solve_alg = _make_solver(yorg, load_pq, n, nload,
                             initial_guess=algebraic_guess,
                             epu=np.asarray(preset.epu, dtype=float))

    # 一致初始化: 在初始发电机角上解代数, 作为流形起点.
    solve_alg.reset()
    z0 = solve_alg(np.asarray(delta0, dtype=float))

    def rhs(t, x):
        dg, om = x[:n], x[n:]
        z = solve_alg(dg)  # 连续法: 约束在每个求值点严格满足
        return _strict_spm_generator_rhs(dg, om, state, preset, basevalue, z)

    # ``tunit`` is a requested output spacing, so include both endpoints:
    # N intervals require N+1 samples.  The previous N-sample grid drifted
    # every checkpoint by one interval over a long fault trajectory.
    steps = max(2, int(round(tlength / tunit)) + 1)
    t_eval = np.linspace(0.0, tlength, steps)
    sol = solve_ivp(rhs, [0.0, tlength], np.r_[delta0, omega0], method=method,
                    t_eval=t_eval, rtol=rtol, atol=atol, dense_output=False)

    delta = sol.y[:n].T
    omega = sol.y[n:].T
    # 沿输出时刻重建代数量(顺序推进 -> 连续法跟踪).
    solve_alg.reset(z0)
    algebraic = np.array([solve_alg(delta[k]) for k in range(sol.t.size)])
    msum = np.sum(preset.m)
    delta_coi = delta - (delta @ preset.m / msum)[:, None]
    return {"time": sol.t, "delta": delta, "omega": omega, "algebraic": algebraic,
            "delta_coi": delta_coi, "success": bool(sol.success), "method": method}


def remap_algebraic_state(z: np.ndarray, from_state, to_state, ngen: int,
                          fallback: np.ndarray | None = None) -> np.ndarray:
    """Map ``[theta_net, voltage_net]`` between bus orderings by bus number.

    This is used at prefault/fault and fault/postfault boundaries.  A removed
    fault bus is omitted from the target state; no zero is inserted as a
    physical network value.  The six-slot MATLAB placeholder is added only by
    reference exporters, never to the Python DAE solver.
    """

    z = np.asarray(z, dtype=float).reshape(-1)
    from_buses = np.asarray(from_state.metadata["transform"], dtype=int)[ngen:]
    to_buses = np.asarray(to_state.metadata["transform"], dtype=int)[ngen:]
    n_from, n_to = from_buses.size, to_buses.size
    if z.size != 2 * n_from:
        raise ValueError("algebraic state width does not match source network")
    theta_from, voltage_from = z[:n_from], z[n_from:]
    fallback_theta = fallback_voltage = None
    if fallback is not None:
        fallback = np.asarray(fallback, dtype=float).reshape(-1)
        if fallback.size != 2 * n_to:
            raise ValueError("fallback algebraic state width does not match target network")
        fallback_theta, fallback_voltage = fallback[:n_to], fallback[n_to:]
    theta_to = np.empty(n_to, dtype=float)
    voltage_to = np.empty(n_to, dtype=float)
    for j, bus in enumerate(to_buses):
        matches = np.flatnonzero(from_buses == int(bus))
        if matches.size == 1:
            i = int(matches[0])
            theta_to[j] = theta_from[i]
            voltage_to[j] = voltage_from[i]
        elif matches.size == 0 and fallback is not None:
            # A faulted bus can re-enter the postfault network.  Carry its
            # physically solved target-network value from ``fallback`` rather
            # than inserting a zero or selecting a cold-start root.
            theta_to[j] = fallback_theta[j]
            voltage_to[j] = fallback_voltage[j]
        else:
            raise ValueError(f"target network bus {int(bus)} missing or duplicated in source")
    return np.r_[theta_to, voltage_to]


def simulate_spm_trajectory(static, *, clear_time: float = 0.2,
                            postfault_time: float = 0.3, tunit: float = 1e-3,
                            method: str = "RK45") -> SpmTrajectoryResult:
    """运行 prefault/fault-on/clearing/postfault recovery 四阶段 SPM 轨迹。

    清除时刻沿用故障段最后一个状态；网络代数状态在每个输出点重新校正并报告残差。
    ``RK45`` 与 ``Radau`` 可用同一接口重复运行比较。
    """

    if clear_time <= 0 or postfault_time <= 0 or tunit <= 0:
        raise ValueError("clear_time, postfault_time and tunit must be positive")
    preset, base = static.preset, static.basevalue
    delta0 = np.asarray(static.prefault.sep_delta, dtype=float)
    omega0 = np.full(preset.ngen, float(static.prefault.sep_omegapu) * base.omega_b)
    # Solve the prefault SEP network first and map its physical branch into
    # the fault network.  Starting fault DAE from [0, ..., 1, ...] can select
    # a different algebraic voltage branch even when both roots have tiny
    # residuals.
    from .spm_energy import solve_spm_network
    pref_state, pref_ok, pref_residual = solve_spm_network(
        delta0, np.asarray(static.prefault.metadata["yfull_mod"], dtype=complex), preset.epu
    )
    if not pref_ok:
        raise RuntimeError(f"prefault SPM network warm-start failed: residual={pref_residual:g}")
    fault_guess = remap_algebraic_state(pref_state, static.prefault, static.fault, preset.ngen)
    fault = simulate_spm_dae(clear_time, tunit, static.fault, preset, base,
                              delta0, omega0, method=method,
                              algebraic_guess=fault_guess)
    delta_clear = fault["delta"][-1]
    omega_clear = fault["omega"][-1]
    post_sep_state, post_sep_ok, post_sep_residual = solve_spm_network(
        np.asarray(static.postfault.sep_delta, dtype=float),
        np.asarray(static.postfault.metadata["yfull_mod"], dtype=complex),
        preset.epu,
    )
    if not post_sep_ok:
        raise RuntimeError(f"postfault SPM network warm-start failed: residual={post_sep_residual:g}")
    post_guess = remap_algebraic_state(
        fault["algebraic"][-1], static.fault, static.postfault, preset.ngen,
        fallback=post_sep_state,
    )
    post = simulate_spm_dae(postfault_time, tunit, static.postfault, preset, base,
                             delta_clear, omega_clear, method=method,
                             algebraic_guess=post_guess)

    # 保留 t=0 的 prefault 标签，故障段从第二个点开始，清除点独立标注。
    f_slice = slice(1, None)
    p_slice = slice(1, None)
    time = np.r_[0.0, np.asarray(fault["time"])[f_slice],
                 clear_time + np.asarray(post["time"])[p_slice]]
    delta = np.vstack([fault["delta"][0], fault["delta"][f_slice], post["delta"][p_slice]])
    omega = np.vstack([fault["omega"][0], fault["omega"][f_slice], post["omega"][p_slice]])
    # fault-on 网络可能移除故障母线，代数量维度小于 postfault；统一输出宽度时
    # 用 NaN 标记“该阶段不存在的节点”，不以 0 冒充可验证状态。
    raw_alg = [np.atleast_2d(fault["algebraic"][0]), fault["algebraic"][f_slice], post["algebraic"][p_slice]]
    nload_parts = [int(a.shape[1] // 2) for a in raw_alg]
    max_nload = max(nload_parts)
    max_alg = 2 * max_nload
    algebraic_parts = []
    theta_parts = []
    voltage_parts = []
    for part, nload_part in zip(raw_alg, nload_parts):
        padded = np.full((part.shape[0], max_alg), np.nan, dtype=float)
        padded[:, :nload_part] = part[:, :nload_part]
        padded[:, max_nload:max_nload + nload_part] = part[:, nload_part:]
        algebraic_parts.append(padded)
        theta_pad = np.full((part.shape[0], max_nload), np.nan, dtype=float)
        voltage_pad = np.full((part.shape[0], max_nload), np.nan, dtype=float)
        theta_pad[:, :nload_part] = part[:, :nload_part]
        voltage_pad[:, :nload_part] = part[:, nload_part:]
        theta_parts.append(theta_pad)
        voltage_parts.append(voltage_pad)
    algebraic = np.vstack(algebraic_parts)
    n_fault = max(0, fault["time"].size - 1)
    n_post = max(0, post["time"].size - 1)
    labels = np.asarray(["prefault"] + ["fault-on"] * max(0, n_fault - 1) +
                        (["clearing"] if n_fault else []) +
                        ["postfault recovery"] * n_post, dtype=object)
    from .spm_energy import spm_network_residual
    yorg_f, load_f, n, nload_f = _algebraic_context(static.fault, preset)
    yorg_p, load_p, _, nload_p = _algebraic_context(static.postfault, preset)
    residual = []
    for k in range(delta.shape[0]):
        if k <= n_fault:
            z = np.asarray(fault["algebraic"][k], dtype=float)
            value = float(np.linalg.norm(spm_network_residual(
                z, delta[k], yorg_f, np.asarray(preset.epu, dtype=float), load_f)))
            residual.append(value if value < 1e-7 else np.nan)
        else:
            post_index = k - n_fault
            raw_post = post["algebraic"][post_index]
            z = np.asarray(raw_post, dtype=float)
            value = float(np.linalg.norm(spm_network_residual(
                z, delta[k], yorg_p, np.asarray(preset.epu, dtype=float), load_p)))
            residual.append(value if value < 1e-7 else np.nan)
    theta_net = np.vstack(theta_parts)
    voltage_net = np.vstack(voltage_parts)
    return SpmTrajectoryResult(time=time, delta_gen=delta, omega_gen=omega,
                               theta_net=theta_net, voltage_net=voltage_net,
                               algebraic_residual=np.asarray(residual),
                               phase_labels=labels,
                               converged=bool(fault["success"] and post["success"] and
                                              np.all(np.isfinite(residual))))


def select_spm_checkpoints(result: SpmTrajectoryResult, *, clear_time: float,
                           postfault_time: float, tunit: float) -> list[dict]:
    """Return fixed, nearest-output SPM checkpoints for cross-validation.

    No interpolation or zero padding is performed.  If a requested time is
    outside the simulated horizon, the record is marked ``available=False``
    instead of silently reusing the final state.  The pre-clearing checkpoint
    is the last output strictly before ``clear_time``; the clearing checkpoint
    must be exactly at the fault trajectory boundary for a valid comparison.
    """
    if clear_time <= 0 or postfault_time <= 0 or tunit <= 0:
        raise ValueError("clear_time, postfault_time and tunit must be positive")
    requested = [
        ("t0", 0.0),
        ("pre-clearing", max(0.0, clear_time - tunit)),
        ("clearing", clear_time),
        ("post-clearing-10ms", clear_time + 0.01),
        ("post-clearing-50ms", clear_time + 0.05),
        ("post-clearing-100ms", clear_time + 0.10),
        ("final", clear_time + postfault_time),
    ]
    times = np.asarray(result.time, dtype=float)
    out: list[dict] = []
    for label, target in requested:
        if times.size == 0 or target < times[0] - 1e-12 or target > times[-1] + 1e-12:
            out.append({"label": label, "requested_time": float(target), "available": False})
            continue
        idx = int(np.argmin(np.abs(times - target)))
        actual = float(times[idx])
        record = {
            "label": label,
            "requested_time": float(target),
            "actual_time": actual,
            "index": idx,
            "available": True,
            "phase": str(result.phase_labels[idx]),
            "delta_gen": np.asarray(result.delta_gen[idx], dtype=float).copy(),
            "omega_gen": np.asarray(result.omega_gen[idx], dtype=float).copy(),
            "theta_net": np.asarray(result.theta_net[idx], dtype=float).copy(),
            "voltage_net": np.asarray(result.voltage_net[idx], dtype=float).copy(),
            "algebraic_residual": float(result.algebraic_residual[idx]),
        }
        out.append(record)
    return out
