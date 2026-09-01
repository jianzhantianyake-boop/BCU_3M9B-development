"""脚本用途:生成暂态轨迹与能量曲线并画图(可运行示例).

使用方法:
    在 ``python_bcu/`` 目录下执行 ``python draw_plots.py``;脚本会做静态初始化,积分
    一段故障轨迹,计算能量,然后弹出两张图并同时保存到 ``figures/`` 目录.
"""

from pathlib import Path
import sys

# 步骤1:把包根目录加入 sys.path,免安装即可导入 bcu_3m9b.
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import numpy as np

from bcu_3m9b import build_static_result
from bcu_3m9b.dynamics import integrate_reduced
from bcu_3m9b.energy import trajectory_energy
from bcu_3m9b.plotting import plot_trajectory, plot_energy


def main() -> None:
    """算数据 -> 画图 -> 存图.

    步骤:
        (1) 静态初始化;(2) 以预故障 SEP 为初值积分一段故障轨迹;(3) 计算能量;
        (4) 分别画轨迹图和能量图,保存到 figures/ 并显示.
    """

    # 步骤1:静态初始化,拿到网络工况与 SEP.
    static = build_static_result()

    # 步骤2:以预故障 SEP 为初值,积分 0.5 s 故障轨迹.
    delta0 = static.prefault.sep_delta
    omega0 = np.full(static.preset.ngen,
                     static.prefault.sep_omegapu * static.basevalue.omega_b)
    traj = integrate_reduced(0.5, 1e-3, static.fault, static.preset,
                             static.basevalue, delta0, omega0)

    # 步骤3:计算这条轨迹上的势能/动能/总能量.
    energy = trajectory_energy(traj, static.preset, static.postfault)

    # 步骤4:画图并保存到 figures/(同时弹窗显示).
    outdir = ROOT / "figures"
    outdir.mkdir(exist_ok=True)
    plot_trajectory(traj, show=False, save=str(outdir / "trajectory.png"))
    plot_energy(energy, show=False, save=str(outdir / "energy.png"))
    print(f"Saved figures to: {outdir}")

    # 如需弹窗查看,取消下面两行注释(需要图形界面).
    # import matplotlib.pyplot as plt
    # plt.show()


if __name__ == "__main__":
    main()
