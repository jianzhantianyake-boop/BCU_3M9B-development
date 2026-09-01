# Python 学习与运行指南

## 1. 环境

核心要求：Python 3.10 或更高版本、NumPy 1.24 或更高版本。当前工作区验证
使用的是 Codex bundled Python 3.12 和 NumPy；当前环境没有 SciPy/Matplotlib。
因此 `run_static.py`、`run_experiment.py`、`run_spm.py`、`smoke_test.py` 都
不应要求 SciPy 或绘图库。

## 2. 第一条命令

```powershell
cd "C:\Users\WangS\Desktop\26年第二学期学习文件\论文文献\literature_ntu_chatgpt_only\代码\BCU_3M9B-main\python_bcu"
& "C:\Users\WangS\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" scripts\run_static.py
```

成功标志是打印 9 个母线、3 台发电机、3 个负荷，三个网络的 `Yred` 维度均
为 `(3, 3)`，并给出预故障/故障后 SEP 残差。残差接近 0 表示本地 Python
非线性求解器完成了当前工况的平衡点计算。

## 3. 逐层学习

### 两机模型

```powershell
& "C:\Users\WangS\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" scripts\run_two_machine.py
```

先读 `bcu_3m9b/two_machine.py` 中的 `electrical_power`、`f_2m`、
`f_2m_fault`、`f_2m_reduce`。状态分别是功角差、相对速度、共同速度；
所有参数通过 `TwoMachineParameters` 传入，不再通过全局变量隐藏依赖。

### 9-bus 静态初始化

```powershell
& "C:\Users\WangS\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" scripts\run_static.py
```

建议依次查看 `cases.py`、`powerflow.py`、`network.py`、`equilibrium.py` 和
`bcu.py`。对应 MATLAB 学习链是 `case9_v2 -> Fun_ResultBack ->
Fun_RXB2Yfull -> Fun_Yfull2Yred -> Fun_SEPfslove`。

### Y-bus 与网络约简

`rxb_to_yfull` 根据线路 `r/x/b` 构造完整导纳矩阵；`add_load_admittance`
把潮流负荷按 `Y=P/V^2-jQ/V^2` 转为恒阻抗；`kron_reduce` 按发电机节点在
前的顺序生成 `Ynn/Ynr/Yrn/Yrr`，并执行
`Yred=Ynn-Ynr*solve(Yrr,Yrn)`。

### SEP、MGP、CUEP

`solve_sep` 的未知量是前 `ngen-1` 台机相对最后一台机的角度和共同速度；
结果再转换为 m 加权 COI 坐标。`find_mgp` 返回 MGP 候选点、每一步残差和
是否按局部极小范数判据找到。`solve_cuep_from_guess` 用 MGP 作为牛顿初值，
并返回残差和是否收敛；若返回 SEP 附近，平台会明确标记为“MGP 候选”，不会
把它强行标成严格 CUEP。

### 能量与 CCT

`potential_energy` 返回线性机械项、无损网络项和耗散路径项；`trajectory_energy`
计算动能/势能/总能量；`energy_cct` 查找总能量越过临界能量的离散步。时域
网格 CCT 由 `time_domain_cct` 完成，稳定性沿用最大相角差 `<2pi` 和末端
接近故障后 SEP 的教学判据。

### 结构保持模型

```powershell
& "C:\Users\WangS\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" scripts\run_spm.py
```

`spm.py` 不使用 SciPy 的 DAE 求解器；每个时间步先求负荷节点角度和电压，
再用发电机摆动方程推进。它用于学习 DAE 的变量分层和网络节点物理意义，
不是对 MATLAB `ode15s` 的严格逐点复刻。

## 4. 结果检查

```powershell
& "C:\Users\WangS\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" tests\smoke_test.py
```

冒烟测试只检查：39-bus 数据形状、9-bus 静态初始化、Yred 维度、SEP 残差和
SPM 代数方程。它不等于 MATLAB/Python 交叉验证；交叉验证需要另存 MATLAB
基准结果、统一参数、统一时间步和误差判据，留到后续阶段。

## 5. 常见问题

* `ModuleNotFoundError: bcu_3m9b`：请从 `python_bcu` 目录运行脚本，或使用
  脚本自带的绝对路径加载逻辑。
* `Python AC 潮流未收敛`：先确认使用项目默认 `case9_v2`，再减小负荷变化，
  或检查 `powerflow.py` 的初值和母线类型。
* `SPM 负荷代数方程未收敛`：先缩短仿真时间，检查 `Preset.s_load`、故障母线
  和网络切线是否导致电压崩溃。
* 没有 Matplotlib：数值计算仍可运行；只有调用 `plotting.py` 才需要额外安装。
