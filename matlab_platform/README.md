# BCU_3M9B：BCU / Energy Function / Transient Stability 教学代码

本项目包含 9-bus、3-generator 系统上的 BCU（Boundary Controlling Unstable
equilibrium point）相关 MATLAB 研究实现，包括网络约简模型、结构保持模型、
能量函数、CCT 和数值暂态仿真。

## 从这里开始

请先阅读 [BCU 学习与运行指南](docs/BCU学习与运行指南.md)。它说明了路径初始化、
推荐阅读顺序、主脚本依赖关系、变量单位、预期输出和常见故障。

MATLAB 中推荐先执行：

```matlab
paths = setup_bcu_paths();
```

旧教程中的 `Y.m` 仍然保留作为兼容入口，但新代码应优先使用
`setup_bcu_paths()`，因为它会根据当前仓库位置自动定位本地 MATPOWER。

## 主要入口

网络约简模型：

- `B3_MM/Cal_MM_CCT.m`：BCU、MGP、CUEP 和能量 CCT；
- `B3_MM/NumSim_MM_Gridframe.m`：三阶段数值暂态仿真与能量后处理；
- `B3_MM/Statable_Region.m`：二维稳定区域探索。

结构保持模型：

- `B3_MM/Cal_MM_CCT_SPM.m`：结构保持模型的 BCU/CCT；
- `B3_MM/NumSim_MM_Gridframe_SPM.m`：结构保持暂态仿真；
- `B3_MM/Statable_Region_SPM.m`：结构保持稳定区域探索。

## 代码边界

`C1_Matpower/matpower7.1` 中的 MATPOWER 官方源码保持原样。项目自有接口、
自定义 case、`Y.m` 和 `B3_MM` 文件增加了中文教学说明；具体修改见
[修改日志](docs/修改日志.md)。

MATLAB 教学化和运行基线见上文；Python 实验仿真平台见下文。当前仍未进行
MATLAB/Python 定量交叉验证。

## Python 实验仿真平台

Python 版本位于 [python_bcu](python_bcu/README_CN.md)。它使用 NumPy 显式传递
案例、潮流、网络约简、SEP、轨迹和能量数据，不依赖 MATLAB 的 base workspace。
当前最小环境只需要 Python 3.10+ 和 NumPy；Matplotlib 是可选绘图库，SciPy
不是运行核心链路的前置条件。

在 PowerShell 中运行：

```powershell
cd "C:\Users\WangS\Desktop\26年第二学期学习文件\论文文献\literature_ntu_chatgpt_only\代码\BCU_3M9B-main\python_bcu"
& "C:\Users\WangS\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" scripts\run_static.py
& "C:\Users\WangS\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" scripts\run_experiment.py
& "C:\Users\WangS\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" tests\smoke_test.py
```

Python 平台当前覆盖：9-bus/39-bus 案例数据、纯 NumPy AC 潮流、恒阻抗负荷、
Kron 约简、两机模型、3M9B 约简模型、SEP/MGP/CUEP 教学接口、势能/动能、
固定步长暂态仿真、能量 CCT、网格 CCT 和结构保持模型的可运行近似。严格
MATLAB/Python 逐点交叉验证留到下一阶段，不把当前冒烟测试称为交叉验证。

详细命令和模块边界见 [Python 学习与运行指南](python_bcu/docs/PYTHON学习与运行指南_CN.md)
以及 [MATLAB-Python 模块映射](python_bcu/docs/MATLAB_PYTHON模块映射_CN.md)。
