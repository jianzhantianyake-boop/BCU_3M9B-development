"""脚本用途：运行 9 母线静态初始化并打印学习所需的关键量。

使用方法：
    命令行执行 ``python scripts/run_static.py``，无需参数；脚本会构造静态结果
    （潮流 + 预故障/故障/故障后网络 + SEP）并打印结果。
"""

from pathlib import Path
import sys

# 步骤1：把包根目录加入 sys.path，免安装即可导入 ``bcu_3m9b``。
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from bcu_3m9b import build_static_result


def main() -> None:
    """执行对应 MATLAB ``Cal_MM_Static`` 的 Python 入口。

    步骤：
        (1) 构造静态结果；(2) 依次打印案例规模、母线电压、预/故障后 SEP 的角度、
        速度和残差，以及各约简导纳维度。
    """

    # 步骤1：一次性构造潮流、网络工况和 SEP。
    result = build_static_result()
    # 步骤2：打印学习所需的各项量（输出全英文）。
    print("=== BCU Python static initialization ===")
    print(f"Case: {result.case.name}")
    print(f"Buses={result.pfdata.nbus}, generators={result.pfdata.ngen}, loads={result.pfdata.nload}")
    print(f"Power-flow residual is bounded inside the static solver; generator buses={result.pfdata.gen_no.tolist()}")
    print("Bus voltage magnitude [pu]:", np.array2string(result.pfdata.voltage[:, 0], precision=6))
    print("Bus voltage angle [deg]:", np.array2string(result.pfdata.voltage[:, 1], precision=6))
    print("Pre-fault SEP angle [rad]:", np.array2string(result.prefault.sep_delta, precision=9))
    print("Post-fault SEP angle [rad]:", np.array2string(result.postfault.sep_delta, precision=9))
    print("Pre-fault SEP speed [pu]:", f"{result.prefault.sep_omegapu:.12g}")
    print("Post-fault SEP speed [pu]:", f"{result.postfault.sep_omegapu:.12g}")
    print("Pre-fault SEP residual:", np.array2string(result.prefault.sep_perr, precision=3))
    print("Post-fault SEP residual:", np.array2string(result.postfault.sep_perr, precision=3))
    print(f"Yred dimensions: pre={result.prefault.yred.shape}, fault={result.fault.yred.shape}, post={result.postfault.yred.shape}")
    print("Note: the Python version does not yet claim point-by-point cross-validation with MATLAB; the current goal is a runnable core pipeline.")


if __name__ == "__main__":
    main()
