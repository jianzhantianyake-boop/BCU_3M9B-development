"""脚本用途：运行两机模型和局部平衡点扫描。

使用方法：
    命令行执行 ``python scripts/run_two_machine.py``，无需参数；脚本装入两机等值
    参数、扫描平衡点并积分一条轨迹。
"""

from pathlib import Path
import sys

# 步骤1：把包根目录加入 sys.path，免安装即可导入。
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from bcu_3m9b.two_machine import TwoMachineParameters, equilibria, simulate_two_machine


def main() -> None:
    """执行两机模型入门实验。

    步骤：
        (1) 由等值阻抗算互导纳并装入参数；(2) 扫描平衡点；(3) 积分一条轨迹并打印
        平衡点数、前几个平衡点和末端状态。
    """

    z1, z2, zl = 0.5j, 0.3j, 0.2844 + 0.0306j
    y12 = 1 / (z1 + z2 + z1 * z2 / zl)
    y1 = 1 / (z1 + zl + z1 * zl / z2)
    y2 = 1 / (z2 + zl + z2 * zl / z1)
    p = TwoMachineParameters(y1.real, y2.real, -y12.real, -y12.imag,
                             1.0, 1.0, 2.2, 1.3, 0.6, 0.5, 0.6, 0.45,
                             (1 / z1).real, (1 / z2).real)
    eps = equilibria(p)
    time, states = simulate_two_machine(p, np.array([0.0, 0.0, 0.0]), 1.0, 0.001, 0.2)
    # 打印实验结果（输出全英文）。
    print("=== Two-machine model Python experiment ===")
    print("Number of equilibria found:", len(eps))
    for item in eps[:5]:
        print("Equilibrium:", item["x"], "non-negative eigenvalues:", item["unstable_dimension"])
    print("Trajectory shape:", states.shape, "final state:", states[-1])


if __name__ == "__main__":
    main()
