"""稳定域与向量场实验工具使用说明。

使用方法：
    算同步坐标功率残差范数用 ``vectorfield_norm``；搜索并分类两机平衡点用
    ``classify_two_machine_equilibria``；对一批初值批量积分并标注稳定/发散用
    ``simulate_grid``；绘图交给可选的 plotting 模块。

对应关系：
    对 MATLAB ``Statable_Region*.m`` 和 ``vectorfield_cal.m`` 的平台化整理。
"""

from __future__ import annotations

from typing import Callable, Dict, Tuple

import numpy as np

from .numerics import newton_solve, numerical_jacobian
from .two_machine import TwoMachineParameters, f_2m


def vectorfield_norm(theta: np.ndarray, yred: np.ndarray, preset) -> float:
    """计算同步坐标下的功率残差范数（对应 vectorfield_cal）。

    使用方法：
        传入功角、约简导纳和参数，返回 COI 功率残差的二范数；范数越小越接近平衡点。
    """

    from .equilibrium import electrical_power

    pe = electrical_power(theta, yred, preset.epu)
    pcoi = np.sum(preset.pmpu - pe)
    return float(np.linalg.norm(preset.pmpu - pe - pcoi / np.sum(preset.m) * preset.m))


def classify_two_machine_equilibria(p: TwoMachineParameters,
                                    angle_grid: np.ndarray | None = None) -> list[dict]:
    """搜索两机三状态平衡点并按非负特征值个数分类。

    使用方法：
        传入两机参数和可选角度网格，返回平衡点列表（转调 two_machine.equilibria）。
    """

    return __import__("bcu_3m9b.two_machine", fromlist=["equilibria"]).equilibria(p, angle_grid)


def simulate_grid(rhs: Callable[[np.ndarray], np.ndarray], starts: np.ndarray,
                  tlength: float = 10.0, tunit: float = 1e-2,
                  divergence: float = 2 * np.pi) -> Dict[str, np.ndarray]:
    """对一批初值批量积分并粗略标注稳定/发散。

    使用方法：
        传入右端函数 ``rhs``、初值集合 ``starts``、总时长、步长和发散阈值，返回
        字典 {starts, terminal, unstable}；状态范数超过阈值即标记该初值为不稳定。
    """

    starts = np.asarray(starts, dtype=float)
    steps = max(2, int(round(tlength / tunit)))
    flags = np.zeros(starts.shape[0], dtype=int)
    terminal = np.zeros_like(starts)
    for i, start in enumerate(starts):
        x = start.copy()
        for _ in range(steps - 1):
            k1 = rhs(x)
            k2 = rhs(x + 0.5 * tunit * k1)
            k3 = rhs(x + 0.5 * tunit * k2)
            k4 = rhs(x + tunit * k3)
            x = x + tunit * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            if not np.all(np.isfinite(x)) or np.linalg.norm(x) > divergence:
                flags[i] = 1
                break
        terminal[i] = x
    return {"starts": starts, "terminal": terminal, "unstable": flags}
