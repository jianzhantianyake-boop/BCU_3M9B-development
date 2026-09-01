"""脚本用途：不依赖 pytest 的平台冒烟测试。

使用方法：
    命令行执行 ``python tests/smoke_test.py``；全部断言通过时打印
    ``smoke_test: PASS``，否则在失败处抛出 AssertionError。
"""

from pathlib import Path
import sys

# 步骤1：把包根目录加入 sys.path，免安装即可导入。
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from bcu_3m9b import build_static_result
from bcu_3m9b.cases import case39_modified
from bcu_3m9b.powerflow import solve_power_flow
from bcu_3m9b.spm import solve_algebraic


def main() -> None:
    """检查案例、静态初始化、SEP 残差和 SPM 代数方程。

    步骤：
        (1) 校验 39 母线案例规模并求潮流；(2) 做 9 母线静态初始化并检查规模与 SEP
        残差；(3) 校验 SPM 负荷代数方程收敛；(4) 全部通过则打印 PASS。
    """

    case39 = case39_modified()
    assert case39.bus.shape == (39, 13)
    case39_pf = solve_power_flow(case39, tol=1e-8)
    assert case39_pf.success and case39_pf.residual_norm < 1e-6
    static = build_static_result()
    assert static.pfdata.nbus == 9
    assert static.pfdata.ngen == 3
    assert static.prefault.yred.shape == (3, 3)
    assert np.linalg.norm(static.prefault.sep_perr) < 1e-6
    assert np.linalg.norm(static.postfault.sep_perr) < 1e-6
    _, ok, residual = solve_algebraic(static.postfault.sep_delta, static.postfault, static.preset)
    assert ok and residual < 1e-6
    print("smoke_test: PASS")


if __name__ == "__main__":
    main()
