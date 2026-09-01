# -*- coding: utf-8 -*-
"""六个实验模式的编排与绘图(对应 MATLAB EXPERIMENT_MODE 的单机/SPM 部分).

使用方法:
    六个入口函数各对应一个 MATLAB 模式,均可直接调用,返回结果字典并把图保存到
    ``figures/``:
        mode_reduced_cct       网络约简:初始化 -> CCT / CUEP(含 δ2-δ3 相平面图)
        mode_reduced_numerical 网络约简:初始化 -> CCT / CUEP -> 数值轨迹(多图)
        mode_reduced_region    网络约简:初始化 -> 二维稳定域搜索(EP + 分界线)
        mode_spm_cct           结构保持:初始化 -> 时域 CCT(正向 SPM 近似)
        mode_spm_numerical     结构保持:初始化 -> 数值轨迹(多图)
        mode_spm_region        结构保持:初始化 -> 二维稳定域(正向网格分类近似)

与 MATLAB 的差异(诚实标注):
    reduced_region 忠实移植 f_reducedstate 梯度系统的 EP 搜索 + type-1 UEP 稳定流形
    分界线.SPM 三个模式原版用 14 维质量矩阵 DAE 与反向 ode15s;这里用正向 SPM 积分
    的网格分类/时域搜索做可运行近似,不声称与 ode15s 逐点一致.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import numpy as np

from .bcu import build_static_result, run_bcu_experiment
from .dynamics import integrate_reduced
from .energy import trajectory_energy
from .equilibrium import electrical_power
from .numerics import numerical_jacobian
from .spm import simulate_spm
from .types import Preset, NetworkState, StaticResult


# 图片默认保存目录:python_bcu/figures/
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"


def _ensure_static(static: Optional[StaticResult]) -> StaticResult:
    """内部:没传 static 就现做一次静态初始化."""

    return static if static is not None else build_static_result()


def _get_plt():
    """内部:延迟导入 matplotlib,缺库时给出清晰提示."""

    try:
        import matplotlib
        matplotlib.use("Agg")  # 无界面也能存图;有界面时可去掉这行
        import matplotlib.pyplot as plt
        return plt
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(
            "matplotlib is not available; run `pip install matplotlib` to enable plotting"
        ) from exc


# ======================= 网络约简梯度系统(用于稳定域) =======================

def reduced_gradient(deltac: np.ndarray, postfault: NetworkState,
                     preset: Preset) -> np.ndarray:
    """网络约简模型的二维梯度系统右端(移植 MATLAB f_reducedstate).

    使用方法:
        传入二维 COI 角 ``deltac=[δ2c, δ3c]``,故障后工况和参数,返回二维导数;
        机 1 角由 COI 约束 ``δ1c = -(m2·δ2c + m3·δ3c)/m1`` 导出.其零点即平衡点.
    """

    m = preset.m
    d2, d3 = float(deltac[0]), float(deltac[1])
    d1 = -(m[1] * d2 + m[2] * d3) / m[0]
    deltacc = np.array([d1, d2, d3])
    pe = electrical_power(deltacc, postfault.yred, preset.epu)
    pcoi = np.sum(preset.pmpu) - np.sum(pe)
    mt = np.sum(m)
    return np.array([preset.pmpu[1] - pe[1] - pcoi / mt * m[1],
                     preset.pmpu[2] - pe[2] - pcoi / mt * m[2]])


def _canonicalize(xep: np.ndarray, m: np.ndarray) -> np.ndarray:
    """内部:把 (δ2c, δ3c) 平衡点解映射到规范周期胞(复刻 MATLAB 变换)."""

    add = (m[1] * xep[0] + m[2] * xep[1]) / m[0]
    xep23 = np.mod(np.array([xep[0] + add, xep[1] + add]), 2.0 * np.pi)
    coi = (xep23[0] * m[1] + xep23[1] * m[2]) / np.sum(m)
    return xep23 - coi


def find_reduced_equilibria(postfault: NetworkState, preset: Preset,
                            grid_points: int = 21, tol: float = 1e-2) -> list:
    """在二维网格上搜索梯度系统平衡点并按稳定性分类.

    使用方法:
        传入故障后工况和参数;返回平衡点列表,每项含规范坐标 ``xep``,原始解
        ``raw``,雅可比特征值,非负特征值个数 ``flag`` 和稳定特征向量 ``vstable``.
    步骤:
        网格初值逐个牛顿求根 -> 规范化去重 -> 求雅可比特征值 -> flag=非负实部个数;
        flag==1 为 type-1 UEP(用于画稳定域分界线).
    """

    from .numerics import newton_solve

    axis = np.linspace(-1.0, 1.0, grid_points) * 2.0 * np.pi
    eps: list = []
    for a in axis:
        for b in axis:
            sol, ok, _, res = newton_solve(
                lambda z: reduced_gradient(z, postfault, preset),
                np.array([a, b]), tol=1e-10, max_iter=200)
            if not ok or res > tol:
                continue
            xep = _canonicalize(sol, preset.m)
            if any(np.max(np.minimum(np.abs(xep - e["xep"]),
                                     np.abs(2 * np.pi - np.abs(xep - e["xep"])))) < tol
                   for e in eps):
                continue
            jac = numerical_jacobian(lambda z: reduced_gradient(z, postfault, preset), sol)
            lam, vec = np.linalg.eig(jac)
            flag = int(np.sum(np.real(lam) >= 0.0))
            stable_cols = np.real(lam) < 0.0
            vstable = np.real(vec[:, stable_cols]) if stable_cols.any() else None
            eps.append({"xep": xep, "raw": sol, "eig": lam, "flag": flag,
                        "vstable": vstable})
    return eps


def _tile_positions(m: np.ndarray) -> np.ndarray:
    """内部:返回把平衡点平移到相邻周期胞的位移矩阵(复刻 position)."""

    mt = np.sum(m)
    p1 = np.array([1 - m[1] / mt, -m[1] / mt])
    p2 = np.array([-m[2] / mt, 1 - m[2] / mt])
    return np.column_stack([p1, p2])


def _backward_manifold(start: np.ndarray, postfault: NetworkState, preset: Preset,
                       tunit: float = 0.05, steps: int = 700,
                       bound: float = 2.6 * np.pi) -> np.ndarray:
    """内部:反向积分梯度系统,追踪 type-1 UEP 的稳定流形(一条半支)."""

    x = np.array(start, dtype=float)
    pts = [x.copy()]
    for _ in range(steps):
        # 反向:dx/dt = -grad;固定步 RK4.
        f = lambda z: -reduced_gradient(z, postfault, preset)
        k1 = f(x)
        k2 = f(x + 0.5 * tunit * k1)
        k3 = f(x + 0.5 * tunit * k2)
        k4 = f(x + tunit * k3)
        x = x + tunit * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        if not np.all(np.isfinite(x)) or np.max(np.abs(x)) > bound:
            break
        pts.append(x.copy())
    return np.array(pts)


# =============================== 六个实验模式 ===============================

def mode_reduced_cct(static: Optional[StaticResult] = None,
                     save_dir: Optional[Path] = None) -> Dict[str, object]:
    """[reduced_cct] 网络约简:初始化 -> CCT / CUEP,并画 δ2-δ3 相平面图.

    使用方法:可选传入 static;打印退出点/MGP/CUEP/LEA·REA CCT,保存相平面图.
    """

    save_dir = Path(save_dir or FIG_DIR)
    save_dir.mkdir(exist_ok=True)
    static = _ensure_static(static)
    rea = robust_reduced_cct(static)  # 先用干净 SEP 算(run_bcu_experiment 会改写 postfault.SEP)
    res = run_bcu_experiment(static)
    print("[reduced_cct] LEA CCT(能量法) = {:.4g}s   REA CCT(时域/有界判据) = {:.4g}s".format(
        res["lea"].cct, rea))
    print("[reduced_cct] 退出点 index =", res["exit_index"], "  CUEP 来源:", res["cuep_source"])

    plt = _get_plt()
    traj = res["fault_trajectory"]
    sep = static.prefault.sep_delta
    cuep = res["cuep_delta"]
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(traj.thetac[:, 1], traj.thetac[:, 2], "-", color="0.6", lw=1.5, label="fault trajectory")
    ax.plot(sep[1], sep[2], "k.", ms=12, label="pre-fault SEP")
    ax.plot(traj.thetac[res["exit_index"], 1], traj.thetac[res["exit_index"], 2],
            "rx", ms=10, mew=2, label="exit point")
    ax.plot(cuep[1], cuep[2], "mo", ms=8, mfc="none", label="CUEP (closest-UEP)")
    ax.set_xlabel(r"$\delta_2$ (COI, rad)")
    ax.set_ylabel(r"$\delta_3$ (COI, rad)")
    ax.set_title("reduced_cct: fault trajectory & critical points")
    ax.grid(True); ax.legend()
    path = save_dir / "reduced_cct_phaseplane.png"
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    print("[reduced_cct] 已保存图:", path)
    return {"result": res, "figure": str(path)}


def mode_reduced_numerical(static: Optional[StaticResult] = None,
                           fault_time: float = 0.2, postfault_time: float = 3.0,
                           prefault_time: float = 0.2, tunit: float = 1e-3,
                           save_dir: Optional[Path] = None) -> Dict[str, object]:
    """[reduced_numerical] 网络约简:三段数值轨迹(prefault/fault/postfault)+ 多图.

    使用方法:可选传入 static 与各段时长;生成 2x2 多图(角度/速度/相平面/能量).
    """

    save_dir = Path(save_dir or FIG_DIR)
    save_dir.mkdir(exist_ok=True)
    static = _ensure_static(static)
    base, preset = static.basevalue, static.preset

    # 三段固定步长积分:prefault(在 SEP 附近保持)-> fault -> postfault.
    d0 = static.prefault.sep_delta
    w0 = np.full(preset.ngen, static.prefault.sep_omegapu * base.omega_b)
    seg_pre = integrate_reduced(prefault_time, tunit, static.prefault, preset, base, d0, w0)
    seg_flt = integrate_reduced(fault_time, tunit, static.fault, preset, base,
                                seg_pre.theta[-1], seg_pre.omega[-1])
    seg_pos = integrate_reduced(postfault_time, tunit, static.postfault, preset, base,
                                seg_flt.theta[-1], seg_flt.omega[-1])

    t0 = seg_pre.time
    t1 = seg_pre.time[-1] + seg_flt.time
    t2 = t1[-1] + seg_pos.time
    time = np.concatenate([t0, t1, t2])
    thetac = np.vstack([seg_pre.thetac, seg_flt.thetac, seg_pos.thetac])
    omegac = np.vstack([seg_pre.omegac, seg_flt.omegac, seg_pos.omegac])
    energy_post = trajectory_energy(seg_pos, preset, static.postfault)

    print("[reduced_numerical] 段时长 pre/fault/post = {:.3g}/{:.3g}/{:.3g}s,总步数 {}".format(
        prefault_time, fault_time, postfault_time, time.size))

    plt = _get_plt()
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for i in range(preset.ngen):
        axes[0, 0].plot(time, thetac[:, i], label=f"gen {i+1}")
    axes[0, 0].set_title("COI angles vs time"); axes[0, 0].set_xlabel("t (s)")
    axes[0, 0].set_ylabel(r"$\delta$ (rad)"); axes[0, 0].grid(True); axes[0, 0].legend()

    for i in range(preset.ngen):
        axes[0, 1].plot(time, omegac[:, i], label=f"gen {i+1}")
    axes[0, 1].set_title("COI relative speed vs time"); axes[0, 1].set_xlabel("t (s)")
    axes[0, 1].set_ylabel(r"$\omega_c$ (rad/s)"); axes[0, 1].grid(True); axes[0, 1].legend()

    axes[1, 0].plot(thetac[:, 1], thetac[:, 2], "-", lw=1.2)
    axes[1, 0].plot(static.postfault.sep_delta[1], static.postfault.sep_delta[2], "k.", ms=10)
    axes[1, 0].set_title(r"$\delta_2$-$\delta_3$ phase plane")
    axes[1, 0].set_xlabel(r"$\delta_2$ (rad)"); axes[1, 0].set_ylabel(r"$\delta_3$ (rad)")
    axes[1, 0].grid(True)

    axes[1, 1].plot(seg_pos.time, energy_post["ep"], label="Ep")
    axes[1, 1].plot(seg_pos.time, energy_post["ek"], label="Ek")
    axes[1, 1].plot(seg_pos.time, energy_post["total"], label="Total")
    axes[1, 1].set_title("post-fault energy vs time"); axes[1, 1].set_xlabel("t (s)")
    axes[1, 1].grid(True); axes[1, 1].legend()

    path = save_dir / "reduced_numerical.png"
    fig.suptitle("reduced_numerical: 3-segment numerical simulation")
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    print("[reduced_numerical] 已保存图:", path)
    return {"time": time, "thetac": thetac, "omegac": omegac, "figure": str(path)}


def mode_reduced_region(static: Optional[StaticResult] = None,
                        grid_points: int = 21,
                        save_dir: Optional[Path] = None) -> Dict[str, object]:
    """[reduced_region] 网络约简:二维稳定域搜索(EP 分类 + type-1 UEP 分界线).

    使用方法:可选传入 static 与网格密度 ``grid_points``;打印各平衡点分类,画出
    δ2-δ3 平面上的平衡点散点(按 flag 着色)与稳定域分界线(稳定流形).
    """

    save_dir = Path(save_dir or FIG_DIR)
    save_dir.mkdir(exist_ok=True)
    static = _ensure_static(static)
    preset, post = static.preset, static.postfault

    eps = find_reduced_equilibria(post, preset, grid_points=grid_points)
    print(f"[reduced_region] 找到平衡点 {len(eps)} 个(flag=非负特征值个数,0=稳定SEP,1=type-1 UEP)")
    for k, e in enumerate(eps):
        print(f"  EP{k}: xep={np.array2string(e['xep'], precision=3)}  flag={e['flag']}")

    plt = _get_plt()
    position = _tile_positions(preset.m)
    shifts = [np.array([0.0, 0.0]), np.array([2 * np.pi, 0.0]), np.array([0.0, 2 * np.pi]),
              np.array([2 * np.pi, 2 * np.pi]), np.array([0.0, 4 * np.pi]),
              np.array([2 * np.pi, 4 * np.pi])]
    colors = {0: "blue", 1: "red", 2: "black"}

    fig, ax = plt.subplots(figsize=(7, 7))
    for e in eps:
        for sh in shifts:
            xep = e["xep"] - position @ sh
            ax.scatter(xep[0], xep[1], c=colors.get(e["flag"], "green"), s=25)
            if e["flag"] == 1 and e["vstable"] is not None:
                v = e["vstable"][:, 0]
                for sign in (+1.0, -1.0):
                    curve = _backward_manifold(xep + sign * 1e-2 * v, post, preset)
                    if curve.shape[0] > 2:
                        ax.plot(curve[:, 0], curve[:, 1], "k-", lw=1.2)
    ax.set_xlim(-2 * np.pi, 2 * np.pi); ax.set_ylim(-2 * np.pi, 2 * np.pi)
    ax.set_xlabel(r"$\delta_2$ (rad)"); ax.set_ylabel(r"$\delta_3$ (rad)")
    ax.set_title("reduced_region: equilibria (blue=SEP, red=type-1 UEP) & stability boundary")
    ax.grid(True)
    path = save_dir / "reduced_region.png"
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    print("[reduced_region] 已保存图:", path)
    return {"equilibria": eps, "figure": str(path)}


def _spm_stable(delta_coi: np.ndarray, sep_delta: np.ndarray,
                bound: float = 2.0 * np.pi, tail_tol: float = 0.3) -> bool:
    """内部:由 SPM 末端 COI 角判定是否稳定(有界且回到 SEP 附近)."""

    if not np.all(np.isfinite(delta_coi)) or np.max(np.abs(delta_coi[-1])) > bound:
        return False
    return bool(np.linalg.norm(delta_coi[-1] - sep_delta) <= tail_tol)


def _bounded(delta_coi: np.ndarray, sep: np.ndarray, limit: float = 1.5 * np.pi) -> bool:
    """内部:判断轨迹是否有界(未失步)--全程有限且离 SEP 不超过 limit."""

    dc = np.asarray(delta_coi)
    return bool(np.all(np.isfinite(dc)) and np.max(np.abs(dc - sep)) < limit)


def robust_reduced_cct(static: StaticResult, fault_max: float = 0.4,
                       samples: int = 21, postfault_time: float = 2.0,
                       tunit: float = 1e-3) -> float:
    """用"有界性"判据做网络约简时域 CCT(比核心 trajectory_stable 更适合轻阻尼).

    使用方法:
        传入 static, 在 [0, fault_max] 上等距取切除时刻, 故障段+故障后段用约简模型
        积分, 以"故障后轨迹是否有界(不失步)"判稳, 返回最大稳定切除时间.
    """

    preset, base = static.preset, static.basevalue
    fault, post = static.fault, static.postfault
    d0 = static.prefault.sep_delta
    w0 = np.full(preset.ngen, static.prefault.sep_omegapu * base.omega_b)
    stable = 0.0
    for tc in np.linspace(0.0, fault_max, samples):
        if tc <= 0:
            continue
        fseg = integrate_reduced(tc, tunit, fault, preset, base, d0, w0)
        pseg = integrate_reduced(postfault_time, tunit, post, preset, base,
                                 fseg.theta[-1], fseg.omega[-1])
        if _bounded(pseg.thetac, post.sep_delta):
            stable = float(tc)
        else:
            break
    return stable


def _spm_until_diverge(tlength, tunit, state, preset, base, d0, w0):
    """内部:跑 SPM, 若代数解中途发散就逐步缩短时长, 返回最长能收敛的一段."""

    length = float(tlength)
    for _ in range(6):
        try:
            return simulate_spm(length, tunit, state, preset, base, d0, w0)
        except RuntimeError:
            length *= 0.6
    return simulate_spm(max(4 * tunit, 0.02), tunit, state, preset, base, d0, w0)


def mode_spm_cct(static: Optional[StaticResult] = None,
                 fault_max: float = 0.4, samples: int = 17,
                 postfault_time: float = 1.0, tunit: float = 2e-3,
                 save_dir: Optional[Path] = None) -> Dict[str, object]:
    """[spm_cct] 结构保持:时域 CCT 搜索(故障段用约简模型, 故障后段用 SPM).

    使用方法:可选传入 static;在 [0, fault_max] 上等距取切除时刻, 故障段用约简摆动
    模型积分到切除(数值稳健), 再用故障后 SPM 积分并判稳, 返回最大稳定切除时间.
    附参考的能量法 LEA CCT. 注:故障网络删了连接母线, 其 SPM 代数方程病态, 故障段
    改用约简模型是刻意的稳健选择.
    """

    save_dir = Path(save_dir or FIG_DIR)
    save_dir.mkdir(exist_ok=True)
    static = _ensure_static(static)
    preset, base = static.preset, static.basevalue
    fault, post = static.fault, static.postfault

    d0 = static.prefault.sep_delta
    w0 = np.full(preset.ngen, static.prefault.sep_omegapu * base.omega_b)
    # 预先求故障后 SEP 的负荷代数解, 作为每次 SPM 首步热启动初值(避免首步不收敛的假失稳).
    from .spm import solve_algebraic
    sep_guess, _, _ = solve_algebraic(post.sep_delta, post, preset)
    clear_times = np.linspace(0.0, fault_max, samples)
    stable_time = 0.0
    for tc in clear_times:
        if tc <= 0:
            continue
        # 故障段用约简摆动模型积分(数值稳健), 故障后段优先用 SPM 判稳.
        fseg = integrate_reduced(tc, min(tunit, 1e-3), fault, preset, base, d0, w0)
        try:
            pseg = simulate_spm(postfault_time, tunit, post, preset, base,
                                fseg.theta[-1], fseg.omega[-1], guess=sep_guess)
            stable = _bounded(pseg["delta_coi"], post.sep_delta)
        except Exception:  # noqa: BLE001  SPM 代数解偶发数值失败: 回退到约简模型有界性判据
            rseg = integrate_reduced(postfault_time, min(tunit, 1e-3), post, preset, base,
                                     fseg.theta[-1], fseg.omega[-1])
            stable = _bounded(rseg.thetac, post.sep_delta)
        if stable:
            stable_time = float(tc)
        else:
            break

    ref = run_bcu_experiment(static)
    print(f"[spm_cct] SPM 时域 CCT ≈ {stable_time:.4g}s "
          f"(参考:约简能量法 LEA CCT = {ref['lea'].cct:.4g}s)")
    print("[spm_cct] 注:SPM 代数解偶发数值失败处已回退用约简模型有界性判据.")
    return {"spm_cct": stable_time, "reference_lea_cct": ref["lea"].cct,
            "tested_clear_times": clear_times}


def mode_spm_numerical(static: Optional[StaticResult] = None,
                       fault_time: float = 0.1, postfault_time: float = 0.6,
                       tunit: float = 2e-3, save_dir: Optional[Path] = None) -> Dict[str, object]:
    """[spm_numerical] 结构保持:数值轨迹 + 多图(故障段约简, 故障后段 SPM).

    使用方法:可选传入 static;故障段用约简模型积分(数值稳健), 故障后段用 SPM, 生成
    发电机 COI 角, 负荷母线电压, δ2-δ3 相平面共 3 图.
    """

    save_dir = Path(save_dir or FIG_DIR)
    save_dir.mkdir(exist_ok=True)
    static = _ensure_static(static)
    preset, base = static.preset, static.basevalue

    d0 = static.prefault.sep_delta
    w0 = np.full(preset.ngen, static.prefault.sep_omegapu * base.omega_b)
    # 故障段用约简摆动模型(数值稳健), 故障后段用结构保持模型 SPM.
    fseg = integrate_reduced(fault_time, min(tunit, 1e-3), static.fault, preset, base, d0, w0)
    pseg = _spm_until_diverge(postfault_time, tunit, static.postfault, preset, base,
                              fseg.theta[-1], fseg.omega[-1])
    time = np.concatenate([fseg.time, fseg.time[-1] + pseg["time"]])
    delta_coi = np.vstack([fseg.thetac, pseg["delta_coi"]])
    # 负荷母线电压来自故障后 SPM 代数解.
    volt = pseg["algebraic"][:, pseg["algebraic"].shape[1] // 2:]
    volt_time = fseg.time[-1] + pseg["time"]

    print(f"[spm_numerical] 故障段(约简) {fseg.time.size} 步 + 故障后段(SPM) {pseg['time'].size} 步")

    plt = _get_plt()
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for i in range(preset.ngen):
        axes[0].plot(time, delta_coi[:, i], label=f"gen {i+1}")
    axes[0].set_title("SPM generator COI angles"); axes[0].set_xlabel("t (s)")
    axes[0].set_ylabel(r"$\delta$ (rad)"); axes[0].grid(True); axes[0].legend()

    for j in range(volt.shape[1]):
        axes[1].plot(volt_time, volt[:, j], lw=1)
    axes[1].set_title("SPM load-bus voltages"); axes[1].set_xlabel("t (s)")
    axes[1].set_ylabel("V (pu)"); axes[1].grid(True)

    axes[2].plot(delta_coi[:, 1], delta_coi[:, 2], "-", lw=1.2)
    axes[2].plot(static.postfault.sep_delta[1], static.postfault.sep_delta[2], "k.", ms=10)
    axes[2].set_title(r"$\delta_2$-$\delta_3$ phase plane")
    axes[2].set_xlabel(r"$\delta_2$ (rad)"); axes[2].set_ylabel(r"$\delta_3$ (rad)")
    axes[2].grid(True)

    path = save_dir / "spm_numerical.png"
    fig.suptitle("spm_numerical: structure-preserving numerical simulation")
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    print("[spm_numerical] 已保存图:", path)
    return {"time": time, "delta_coi": delta_coi, "figure": str(path)}


def mode_spm_region(static: Optional[StaticResult] = None, grid_points: int = 15,
                    sim_time: float = 0.6, tunit: float = 3e-3,
                    span: float = 0.9 * np.pi,
                    save_dir: Optional[Path] = None) -> Dict[str, object]:
    """[spm_region] 结构保持:二维稳定域(正向网格分类近似 + 解析分界线参考).

    使用方法:可选传入 static 与网格密度;在 δ2-δ3 平面撒初值,用正向 SPM 积分固定
    时长后按"是否有界"判稳,画出稳定/不稳定散点,并叠加约简模型的解析分界线作参考.
    步骤:
        对每个网格初值构造发电机角(机1 由 COI 约束导出)与静息速度 -> SPM 积分 ->
        按有界性分类.这是正向分类近似,非原版反向 DAE 稳定流形.
    """

    save_dir = Path(save_dir or FIG_DIR)
    save_dir.mkdir(exist_ok=True)
    static = _ensure_static(static)
    preset, base, post = static.preset, static.basevalue, static.postfault
    sep = post.sep_delta
    w0 = np.full(preset.ngen, post.sep_omegapu * base.omega_b)
    # SEP 负荷代数解作为每个网格点 SPM 首步的热启动初值.
    from .spm import solve_algebraic
    sep_guess, _, _ = solve_algebraic(sep, post, preset)

    axis = np.linspace(sep[1] - span, sep[1] + span, grid_points)
    axis2 = np.linspace(sep[2] - span, sep[2] + span, grid_points)
    stable_pts, unstable_pts = [], []
    for d2 in axis:
        for d3 in axis2:
            d1 = -(preset.m[1] * d2 + preset.m[2] * d3) / preset.m[0]
            delta0 = np.array([d1, d2, d3])
            # 稳定判据用"有界"而非"回到 SEP":本系统阻尼很轻, 短时不会settle,
            # 稳定点做有界振荡, 不稳定点角度跑飞(或代数解发散).
            try:
                seg = simulate_spm(sim_time, tunit, post, preset, base, delta0, w0,
                                   guess=sep_guess)
                dc = seg["delta_coi"]
                ok = bool(np.all(np.isfinite(dc)) and np.max(np.abs(dc - sep)) < 1.5 * np.pi)
            except Exception:  # noqa: BLE001  代数解发散即失步, 判为不稳定
                ok = False
            (stable_pts if ok else unstable_pts).append((d2, d3))

    print(f"[spm_region] 网格 {grid_points}x{grid_points}:稳定 {len(stable_pts)} 点,"
          f"不稳定 {len(unstable_pts)} 点")

    plt = _get_plt()
    fig, ax = plt.subplots(figsize=(7, 7))
    if stable_pts:
        s = np.array(stable_pts); ax.scatter(s[:, 0], s[:, 1], c="tab:blue", s=30, label="stable")
    if unstable_pts:
        u = np.array(unstable_pts); ax.scatter(u[:, 0], u[:, 1], c="tab:red", s=30, marker="x", label="unstable")
    # 叠加约简模型的解析分界线(type-1 UEP 稳定流形)作为参考.
    try:
        for ep in find_reduced_equilibria(post, preset, grid_points=15):
            if ep["flag"] == 1 and ep["vstable"] is not None:
                v = ep["vstable"][:, 0]
                for sign in (+1.0, -1.0):
                    curve = _backward_manifold(ep["xep"] + sign * 1e-2 * v, post, preset)
                    if curve.shape[0] > 2:
                        ax.plot(curve[:, 0], curve[:, 1], "g-", lw=1.2, alpha=0.7)
        ax.plot([], [], "g-", lw=1.2, label="reduced-model boundary")
    except Exception:  # noqa: BLE001
        pass
    ax.plot(sep[1], sep[2], "k*", ms=14, label="post-fault SEP")
    ax.set_xlim(sep[1] - span, sep[1] + span); ax.set_ylim(sep[2] - span, sep[2] + span)
    ax.set_xlabel(r"$\delta_2$ (rad)"); ax.set_ylabel(r"$\delta_3$ (rad)")
    ax.set_title("spm_region: forward-classified points + reduced-model boundary")
    ax.grid(True); ax.legend()
    path = save_dir / "spm_region.png"
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    print("[spm_region] 已保存图:", path)
    return {"stable": stable_pts, "unstable": unstable_pts, "figure": str(path)}


# ======================= 两机模型三维稳定域(模式 7 / 8) =======================
# 说明: MATLAB Statable_Region_3D.m 与 _GFL.m 结构完全一样, 都用同一套两机完整模型
# f_2m(状态 [δ12, ω12, ω_sum]), 只是参数不同(GFL 版阻抗含电阻, 惯量更低). 因此这里
# 一个函数 + 两组参数即可, 不需要额外的 GFL 动力学模型.

def _params_two_machine(gfl: bool):
    """返回两机模型参数(gfl=False 对应 3D 例, True 对应 GFL 例).

    使用方法: 内部辅助; 由阻抗算互导纳后装成 TwoMachineParameters.
    """

    from .two_machine import TwoMachineParameters

    if not gfl:
        z1, z2, zl = 0.5j, 0.3j, 0.2844 + 0.0306j
        # d2=0.5: 对齐 MATLAB Statable_Region_3D.m 的 EP 实际取值(源码 D2 双重赋值 0.45/0.5,
        # 平衡点用 0.5; 用户 2026-09-01 确认以 0.5 为准).
        pm1, pm2, h1, h2, d1, d2 = 2.2, 1.3, 0.6, 0.5, 0.6, 0.5
    else:
        z1, z2, zl = 0.05 + 0.5j, 0.01 + 0.3j, 0.2844 + 0.0306j
        pm1, pm2, h1, h2, d1, d2 = 1.33, 0.6, 1.0, 1.0, 0.6, 0.6
    y12 = 1 / (z1 + z2 + z1 * z2 / zl)
    y1 = 1 / (z1 + zl + z1 * zl / z2)
    y2 = 1 / (z2 + zl + z2 * zl / z1)
    return TwoMachineParameters(y1.real, y2.real, -y12.real, -y12.imag,
                                1.0, 1.0, pm1, pm2, h1, h2, d1, d2,
                                (1 / z1).real, (1 / z2).real)


def _close_2m(a: np.ndarray, b: np.ndarray, tol: float) -> bool:
    """内部: 判断两机平衡点是否重复(δ12 按 2π 周期比较)."""

    d = np.abs(np.asarray(a) - np.asarray(b))
    d[0] = min(d[0], abs(2 * np.pi - d[0]))
    return bool(np.max(d) < tol)


def _stable_real_basis(lam: np.ndarray, vec: np.ndarray) -> list:
    """内部: 从雅可比特征分解取稳定子空间的实正交基.

    使用方法: 传入特征值 lam 与特征向量矩阵 vec; 返回长度=稳定特征值个数的实向量
    列表(复共轭对用实部/虚部张成实基, 再 QR 正交化).
    """

    stable = np.where(lam.real < 0.0)[0]
    if stable.size == 0:
        return []
    cols, used = [], set()
    for i in stable:
        if i in used:
            continue
        v = vec[:, i]
        if abs(lam[i].imag) < 1e-9:
            cols.append(v.real)
        else:
            cols.append(v.real)
            cols.append(v.imag)
            for j in stable:  # 标记共轭伙伴, 避免重复
                if j != i and abs(lam[j] - np.conj(lam[i])) < 1e-6:
                    used.add(j)
        used.add(i)
    B = np.array(cols[:stable.size], dtype=float).T
    Q, _ = np.linalg.qr(B)
    return [Q[:, k] for k in range(B.shape[1])]


def find_two_machine_equilibria_3d(p, grid_points: int = 11, tol: float = 1e-2) -> list:
    """搜索两机完整模型 f_2m 的平衡点并分类.

    使用方法: 传入参数 p; 以 [δ12_grid, 0, 0] 为初值逐个牛顿求根, 按 δ12 的 2π 周期
    去重, 返回列表, 每项含解 ``xep``, 非负特征值个数 ``flag``, 稳定子空间实基 ``basis``.
    """

    from .two_machine import f_2m
    from .numerics import newton_solve

    eps = []
    for d in np.linspace(0.0, 1.0, grid_points) * 2.0 * np.pi:
        sol, ok, _, res = newton_solve(lambda z: f_2m(z, p),
                                       np.array([d, 0.0, 0.0]), tol=1e-10, max_iter=300)
        if not ok or res > tol:
            continue
        if any(_close_2m(sol, e["xep"], tol) for e in eps):
            continue
        jac = numerical_jacobian(lambda z: f_2m(z, p), sol)
        lam, vec = np.linalg.eig(jac)
        eps.append({"xep": sol, "flag": int(np.sum(lam.real >= 0.0)),
                    "basis": _stable_real_basis(lam, vec), "eig": lam})
    return eps


def _backward_2m(start: np.ndarray, p, steps: int = 300, tunit: float = 0.02,
                 bound=(4.0 * np.pi, 30.0, 220.0)) -> np.ndarray:
    """内部: 反向积分 f_2m(dx/dt = -f_2m), 追踪稳定流形一条曲线, 越界即停."""

    from .two_machine import f_2m

    x = np.array(start, dtype=float)
    pts = [x.copy()]
    for _ in range(steps):
        f = lambda z: -f_2m(z, p)
        k1 = f(x); k2 = f(x + 0.5 * tunit * k1)
        k3 = f(x + 0.5 * tunit * k2); k4 = f(x + tunit * k3)
        x = x + tunit * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        if (not np.all(np.isfinite(x)) or abs(x[0]) > bound[0]
                or abs(x[1]) > bound[1] or abs(x[2]) > bound[2]):
            break
        pts.append(x.copy())
    return np.array(pts)


def mode_two_machine_region_3d(static: Optional[StaticResult] = None, gfl: bool = False,
                               grid_points: int = 11, alpha_count: int = 72,
                               back_steps: int = 300, back_tunit: float = 0.02,
                               perturb: float = 0.5,
                               save_dir: Optional[Path] = None) -> Dict[str, object]:
    """[two_machine_region_3d] 两机完整模型的三维稳定域与稳定流形.

    使用方法:
        独立于电网 static(两机模型自成一体, static 参数被忽略). gfl=True 换成 GFL 参数.
        找平衡点 -> 3D 散点(按 flag 着色) -> 对 type-1 UEP 扫其二维稳定子空间反向积分出
        稳定流形曲面, type-2 UEP 出一维曲线 -> 叠加一条"故障->切除"样例轨迹.
    步骤:
        (1) 以 [δ12,0,0] 网格找 EP 并按非负特征值个数分类; (2) δ12 方向 ±2π 平铺;
        (3) type-1: 沿二维稳定基画 alpha_count 条反向积分曲线(近似流形曲面);
        (4) type-2: 沿一维稳定基正负两条; (5) 画故障+故障后样例轨迹.
    """

    save_dir = Path(save_dir or FIG_DIR)
    save_dir.mkdir(exist_ok=True)
    from .two_machine import simulate_two_machine

    p = _params_two_machine(gfl)
    eps = find_two_machine_equilibria_3d(p, grid_points=grid_points)
    tag = "gfl" if gfl else "3d"
    print(f"[two_machine_region_{tag}] 找到平衡点 {len(eps)} 个 "
          f"(flag: 0=稳定, 1=type-1 UEP, 2=type-2 UEP)")
    for k, e in enumerate(eps):
        print(f"  EP{k}: xep={np.array2string(e['xep'], precision=3)}  flag={e['flag']}")

    plt = _get_plt()
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(projection="3d")
    colors = {0: "blue", 1: "red", 2: "magenta", 3: "green"}
    shifts = [-2.0 * np.pi, 0.0, 2.0 * np.pi]
    for e in eps:
        for sh in shifts:
            xep = e["xep"] + np.array([sh, 0.0, 0.0])
            ax.scatter(xep[0], xep[1], xep[2], c=colors.get(e["flag"], "black"), s=25)
            if e["flag"] == 1 and len(e["basis"]) == 2:
                for a in np.linspace(0.0, 2.0 * np.pi, alpha_count, endpoint=False):
                    vp = e["basis"][0] * np.sin(a) + e["basis"][1] * np.cos(a)
                    curve = _backward_2m(xep + perturb * vp, p, back_steps, back_tunit)
                    if curve.shape[0] > 2:
                        ax.plot(curve[:, 0], curve[:, 1], curve[:, 2], color="0.6", lw=0.4)
            elif e["flag"] == 2 and len(e["basis"]) >= 1:
                for beta in (+1.0, -1.0):
                    curve = _backward_2m(xep + perturb * beta * e["basis"][0], p,
                                         back_steps, back_tunit)
                    if curve.shape[0] > 2:
                        ax.plot(curve[:, 0], curve[:, 1], curve[:, 2], color="0.3", lw=0.7)

    # 故障 -> 切除样例轨迹(从第一个平衡点出发).
    t_fault = 15.0 if gfl else 1.65
    if eps:
        tstep = 2e-3
        time, states = simulate_two_machine(p, eps[0]["xep"], tlength=t_fault + 20.0,
                                            tunit=tstep, fault_until=t_fault)
        nf = int(round(t_fault / tstep))
        ax.plot(states[:nf, 0], states[:nf, 1], states[:nf, 2], "r-", lw=1.5, label="fault")
        ax.plot(states[nf:, 0], states[nf:, 1], states[nf:, 2], "b-", lw=1.2, label="post-fault")

    ax.set_xlim(-2 * np.pi, 2 * np.pi); ax.set_ylim(-15, 10); ax.set_zlim(-100, 100)
    ax.set_xlabel(r"$\delta_{12}$"); ax.set_ylabel(r"$\omega_{12}$"); ax.set_zlabel(r"$\omega_{sum}$")
    ax.set_title(f"two_machine_region_{tag}: 3D stability region & stable manifolds")
    ax.view_init(elev=20, azim=-60)
    ax.legend(loc="upper left")
    path = save_dir / f"two_machine_region_{tag}.png"
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    print(f"[two_machine_region_{tag}] 已保存图: {path}")
    return {"equilibria": eps, "figure": str(path)}


def mode_two_machine_region_3d_gfl(static: Optional[StaticResult] = None,
                                   save_dir: Optional[Path] = None,
                                   **kwargs) -> Dict[str, object]:
    """[two_machine_region_3d_gfl] 同上, 换用 GFL 参数(阻抗含电阻, 惯量更低)."""

    return mode_two_machine_region_3d(static, gfl=True, save_dir=save_dir, **kwargs)


# 模式名 -> 函数,供 main.py 调用.
MODES = {
    "reduced_cct": mode_reduced_cct,
    "reduced_numerical": mode_reduced_numerical,
    "reduced_region": mode_reduced_region,
    "spm_cct": mode_spm_cct,
    "spm_numerical": mode_spm_numerical,
    "spm_region": mode_spm_region,
    "two_machine_region_3d": mode_two_machine_region_3d,
    "two_machine_region_3d_gfl": mode_two_machine_region_3d_gfl,
}
