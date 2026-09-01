# MATLAB-Python 模块映射

Python 平台按“研究核心”和“入口脚本”组织，不机械复制 MATLAB 的 base
workspace 结构。下表给出当前项目自有 MATLAB 文件的主要去向。

| MATLAB 文件/类别 | Python 去向 | 当前状态 |
|---|---|---|
| `case9_v2.m` | `bcu_3m9b/cases.py::case9_v2` | 已转写并运行 |
| `case39_modified.m` | `bcu_3m9b/cases.py::case39_modified` | 数据已转写，默认未作为最小入口 |
| `Fun_ResultBack.m` | `powerflow.py::solve_power_flow/to_pfdata` | 已转写并运行 |
| `Fun_RXB2Yfull.m` | `network.py::rxb_to_yfull` | 已转写 |
| `Fun_Yfull2Yred.m` | `network.py::kron_reduce` | 已转写 |
| `Fun_Yfull2Yfull.m` | `network.py::reorder_structure_preserved` | 已转写 |
| `Fun_SEPfslove.m`、`Fun_SEPcheck.m`、`Fun_SEPiteration.m` | `equilibrium.py`、`matlab_compat.py` | 已转写接口 |
| `f_2m.m`、`f_2m_fault.m`、`f_2m_reduce.m` | `two_machine.py` | 已转写并运行 |
| `F_3M9B_MR_ODE.m` | `dynamics.py::reduced_rhs` | 已转写 |
| `Fun_TrajIter_SRF.m`、`Fun_TrajIter_StableCheck_SRF.m` | `dynamics.py`、`matlab_compat.py` | 已转写 |
| `Fun_Cal_Exitpoint.m`、`Fun_Cal_CCT_Real.m` | `dynamics.py::find_exitpoint/time_domain_cct` | 已转写为批处理接口 |
| `Fun_Cal_PotentialEnergy.m`、`Fun_Cal_CCT_Energy.m` | `energy.py` | 已转写 |
| `Fun_Cal_MGP*.m`、`Fun_Cal_UpdateStartPoint.m` | `energy.py::find_mgp` | 已转写为非交互式候选搜索 |
| `Fun_Cal_DampingEnergy.m` | `energy.py` 的轨迹能量接口 | 核心动/势能已覆盖，完整交互式绘图后续扩展 |
| `Fun_Cal_GenEMF.m` | `powerflow.py::generator_internal_emf` | 已转写 |
| `F_3M9B_SP_ODE.m`、`F_3M9B_SP_DAE.m` | `spm.py` | 已转写为纯 NumPy DAE 近似 |
| `Fun_AEfslove_SPM.m`、`Fun_AEiteration_SPM.m` | `spm.py::algebraic_residual/solve_algebraic` | 已转写 |
| `f_reducedstate*.m` | `dynamics.py`、`stability_region.py` | 约简核心已覆盖 |
| `Statable_Region*.m`、`vectorfield_cal.m` | `stability_region.py`、`plotting.py` | 网格/向量场核心已覆盖 |
| `Plot_3Dstate.m`、`plottmp.m` | `plotting.py` | 改为可选绘图 API |
| `Cal_MM_Static*.m` | `bcu.py::build_static_result` | 已转写并运行 |
| `Cal_MM_CCT*.m`、`NumSim_MM_Gridframe*.m` | `bcu.py::run_bcu_experiment` 和 scripts | 已形成 Python 入口 |

第三方 `C1_Matpower/matpower7.1` 官方源码没有复制和修改。项目自有 MATLAB
接口和 case 数据被作为 Python 平台的数据/算法层重新组织。由于 MATLAB 原
入口依赖交互式 `figure/input`、`evalin('base')` 和 `ode15s`，Python 版本采用
显式参数、批处理返回值和可选绘图；这属于平台化改写，不能在交叉验证尚未建立
前称为严格逐文件等价。
