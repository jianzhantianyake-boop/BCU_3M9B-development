"""脚本用途：运行结构保持模型（SPM）的最小实验。

使用方法：
    命令行执行 ``python scripts/run_spm.py``，无需参数；脚本先做静态初始化，再解
    代数负荷节点电压并推进少量时间步。
"""

from pathlib import Path
import sys

# 步骤1：把包根目录加入 sys.path，免安装即可导入。
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from bcu_3m9b.bcu import build_static_result
from bcu_3m9b.spm import simulate_spm


def main() -> None:
    """求代数负荷节点电压并推进少量时间步。

    步骤：
        (1) 构造静态结果；(2) 以故障后 SEP 为初值跑一小段 SPM；(3) 打印时间步数、
        代数状态维度和末端功角/电压。
    """

    static = build_static_result()
    omega0 = np.full(static.preset.ngen, static.postfault.sep_omegapu * static.basevalue.omega_b)
    data = simulate_spm(0.05, 0.002, static.postfault, static.preset,
                        static.basevalue, static.postfault.sep_delta, omega0)
    # 打印实验结果（输出全英文）。
    print("=== SPM Python minimal experiment ===")
    print("Number of time steps:", data["time"].size)
    print("Algebraic state dimension:", data["algebraic"].shape[1])
    print("Final COI rotor angle:", data["delta_coi"][-1])
    print("Final network voltage:", data["algebraic"][-1, 6:])


if __name__ == "__main__":
    main()
