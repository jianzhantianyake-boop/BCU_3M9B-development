# Python 复写修改日志

## 2026-08-25

### 新增

* `python_bcu/bcu_3m9b/types.py`：显式数据结构，替代 MATLAB base workspace。
* `cases.py`：项目 `case9_v2` 和 `case39_modified` 数据入口。
* `powerflow.py`：纯 NumPy AC 潮流、分支潮流和 `pfdata` 转换。
* `network.py`：Y-bus、恒阻抗负荷、故障母线/线路和 Kron 约简。
* `equilibrium.py`：SEP 残差、COI 变换、Newton 求解和检查。
* `dynamics.py`：3M9B 约简模型、RK4、退出点、稳定性和时域 CCT 网格。
* `energy.py`：势能、动能、MGP、能量轨迹和 LEA CCT。
* `two_machine.py`：两机完整/故障/约简模型和局部平衡点扫描。
* `spm.py`：结构保持模型的可运行纯 NumPy DAE 近似。
* `stability_region.py`、`plotting.py`：稳定域核心和可选绘图接口。
* `scripts/`、`tests/smoke_test.py`、中文运行文档和模块映射。

### 运行行为说明

* 没有修改 MATLAB 文件，也没有修改 MATPOWER 官方源码。
* Python 潮流使用项目自定义 case 数据和自己的有限差分牛顿求解器。
* Python 使用显式数据对象，不再隐式读取 `evalin('base')`。
* SPM 使用“代数方程 Newton + 发电机显式 RK4”的可运行近似，不声称等同于
  MATLAB `ode15s` 的 DAE 数值轨迹。
* MGP/CUEP 和 CCT 已有可批处理入口，但 MATLAB/Python 定量交叉验证尚未实施。

### 验证结果

* `scripts/run_static.py`：已运行，9-bus/3-generator，三个 `Yred` 均为 3x3，
  预故障和故障后 SEP 残差达到约 `1e-11` 量级。
* `scripts/run_experiment.py`：已运行，故障轨迹、退出点、MGP 候选、LEA 和
  REA 网格 CCT 均能返回。
* `scripts/run_spm.py`：已运行，SPM 代数状态和暂态状态均能返回。
* `scripts/run_two_machine.py`：已运行，搜索到多个平衡点并完成轨迹积分。
* `tests/smoke_test.py`：已通过。

### 尚未解决/后续工作

* 统一 MATLAB 基准结果文件与 Python 输出的定量交叉验证。
* 完善 SPM 的故障切换、ZIP 负荷和更严格 DAE 积分器。
* 完整重做交互式阻尼能量分组、三维稳定域曲面和所有科研绘图样式。
* 39-bus 默认动态参数、故障设置和对应的 10 机研究入口尚未作为默认实验。

### 明确未修改

`C1_Matpower/matpower7.1` 下 MATPOWER 官方源码未修改；Python 只读取并转写
项目自有 `case9_v2.m`、`case39_modified.m` 和项目自有接口所需的数据。
