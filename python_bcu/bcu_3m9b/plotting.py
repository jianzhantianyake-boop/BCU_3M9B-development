"""可选绘图接口使用说明。

使用方法：
    安装 Matplotlib 后调用 ``plot_trajectory`` 画 COI 功角/相对速度，调用
    ``plot_energy`` 画动能/势能/总能量；未安装时导入失败会给出明确提示，不会阻塞
    潮流、SEP、能量和轨迹等核心计算。
"""

from __future__ import annotations

from typing import Optional

import numpy as np


def _plt():
    # 延迟导入 Matplotlib：缺库时抛出清晰提示，核心计算不受影响。
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("Matplotlib is not available in this Python environment; core computation is unaffected, install a plotting library separately") from exc
    return plt


def plot_trajectory(traj, title: str = "Reduced 3-machine model trajectory", show: bool = True,
                    save: Optional[str] = None):
    """画 COI 功角和相对速度轨迹。

    使用方法：
        传入 ``Trajectory`` 对象，可选标题、是否显示和保存路径；返回 (fig, axes)。
    """

    plt = _plt()
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
    axes[0].plot(traj.time, traj.thetac)
    axes[0].set_ylabel("COI angle / rad")
    axes[1].plot(traj.time, traj.omegac)
    axes[1].set_ylabel("COI relative speed / rad/s")
    axes[1].set_xlabel("Time / s")
    fig.suptitle(title)
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=160)
    if show:
        plt.show()
    return fig, axes


def plot_energy(energy: dict, title: str = "Energy function", show: bool = True,
                save: Optional[str] = None):
    """画动能、势能和总能量。

    使用方法：
        传入 ``trajectory_energy`` 返回的能量字典，可选标题、是否显示和保存路径；
        返回 (fig, ax)。
    """

    plt = _plt()
    fig, ax = plt.subplots(figsize=(8, 4))
    n = len(energy["total"])
    t = np.arange(n)
    ax.plot(t, energy["ep"], label="Potential energy")
    ax.plot(t, energy["ek"], label="Kinetic energy")
    ax.plot(t, energy["total"], label="Total energy")
    ax.set_title(title)
    ax.set_xlabel("Discrete step")
    ax.legend()
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=160)
    if show:
        plt.show()
    return fig, ax
