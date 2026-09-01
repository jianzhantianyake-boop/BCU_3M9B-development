"""基于 NumPy 的 BCU/能量函数/暂态稳定教学仿真平台。

使用方法：
    从本包直接导入常用入口：案例用 ``case9_v2``/``case39_modified``，潮流用
    ``solve_power_flow``/``to_pfdata``，静态初始化用 ``build_static_result``，
    默认参数用 ``default_preset``，一键实验用 ``run_bcu_experiment``。

设计说明：
    本包把原 MATLAB 工程的核心研究对象改写成显式 Python 数据结构，避免依赖
    MATLAB 的 base workspace，并尽量保持母线顺序、发电机顺序、标幺值约定和 COI
    坐标定义。
"""

from .types import BaseValue, CaseData, Preset, PowerFlowResult, StaticResult
from .cases import case9_v2, case39_modified
from .powerflow import solve_power_flow, to_pfdata
from .bcu import build_static_result, default_preset, run_bcu_experiment

__all__ = [
    "BaseValue",
    "CaseData",
    "Preset",
    "PowerFlowResult",
    "StaticResult",
    "case9_v2",
    "case39_modified",
    "solve_power_flow",
    "to_pfdata",
    "build_static_result",
    "default_preset",
    "run_bcu_experiment",
]
