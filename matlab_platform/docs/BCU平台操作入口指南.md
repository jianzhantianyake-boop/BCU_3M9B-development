# BCU_3M9B 仿真平台操作入口指南

本指南面向"如何运行实验、如何微调参数、如何批量扫参"。核心理念：**所有可调参数集中在
一个文件 `bcu_config.m`，绝不深入修改核心计算脚本**。

---

## 1. 环境要求

| 项 | 要求 |
|----|------|
| MATLAB | R2021b+（本平台在 **R2024a Update 8** 验证通过） |
| 工具箱 | Optimization Toolbox（`fsolve`） |
| MATPOWER | 7.1（放在 `C1_Matpower/matpower7.1`，作为外部依赖） |
| VSCode（可选） | 安装 `MathWorks.language-matlab` 扩展；`.vscode/settings.json` 已配置 `installPath` |

在 VSCode 里打开 `matlab_platform` 文件夹，任意 `.m` 文件右上角有 ▶ **Run**；`Ctrl+Enter` 逐节运行。

---

## 2. 三步快速开始

```matlab
% 1) （首次或换机器时）初始化路径
setup_bcu_paths();

% 2) 编辑 bcu_config.m，改你要的参数，保存

% 3) 运行主入口
run_bcu
```

`run_bcu` 会：**校验配置 → 打印摘要 → 存参数快照到 results/ → 运行实验链路 → 自检**。

---

## 3. 参数总表（都在 `bcu_config.m` 里改）

| 组 | 字段 | 含义 / 单位 | 约束 |
|----|------|------------|------|
| ① 实验 | `mode` | 实验链路 | 见下方模式表 |
| ② 系统 | `CaseName` | MATPOWER case 名 | 文件须存在 |
| | `f_base` | 系统基频 Hz | > 0（9-bus 用 60） |
| ③ 发电机 | `m` | 惯性 M=2H/ωs (pu·s²/rad) | 长度 = 发电机数 |
| | `damping_ratio` | 阻尼比 d/m（实际 d = m·ratio） | 同上 |
| | `Pm` | 机械功率 (pu) | 同上 |
| | `xd1` | 暂态电抗 xd' (pu) | 同上 |
| | `E` | 内电势幅值 (pu) | 同上 |
| ④ 故障 | `faultline` | 故障支路 [From;To] | 须是 case 中存在的支路 |
| | `faultposition` | 0 取 From，1 取 To | 0 或 1 |
| ⑤ 数值 | `EquCal` | 1=Newton, 2=fsolve | 1 或 2 |
| | `PathEnergyCal` | 0=Ray, N=N段梯形, -1=忽略 | ≥ -1 整数 |
| | `Tfault` | 故障积分时长 (s) | > 预计 CCT |
| | `Tunit` | 积分步长 (s) | 0 < Tunit < Tfault |
| | `TolFun`/`TolX` | fsolve 容差 | — |
| ⑥ 输出 | `run_selfcheck` | 跑后自检 | true/false |
| | `save_snapshot` | 存参数快照 | true/false |

### 实验模式（`cfg.mode`）

| mode | 链路 |
|------|------|
| `reduced_cct` | 网络约简：初始化 → CCT/CUEP |
| `reduced_numerical` | 约简：CCT/CUEP → 数值轨迹 |
| `reduced_region` | 约简：二维稳定区域 |
| `spm_cct` / `spm_numerical` / `spm_region` | 结构保持模型对应链路 |
| `two_machine_region_3d[_gfl]` | 独立两机三维示例 |

> 参数注入当前完整覆盖 `reduced_*` 链路；SPM/两机链路暂用其脚本自带参数（见工作日志"下一步"）。

---

## 4. 三种典型操作

### 4.1 微调单个工况
改 `bcu_config.m` 相应字段 → `run_bcu`。填错维度或不存在的支路，会在计算前被 `[ERR]` 拦下。

### 4.2 切换实验
只改 `cfg.mode`，`run_bcu` 自动跑对应链路。

### 4.3 批量扫参
编辑 `run_bcu_sweep.m` 顶部的"扫描定义"三行，然后运行 `run_bcu_sweep`：

```matlab
sweep.param  = 'damping_ratio';                 % 要扫的字段
sweep.values = {[0.05;0.05;0.05], [0.1;0.1;0.1], [0.2;0.2;0.2], [0.3;0.3;0.3]};
sweep.fixed  = struct('Tfault', 0.8);           % 每组附加的固定覆盖（加速）
```

输出：`results/sweep_<时间戳>.{csv,mat}` + 一张 CCT-参数曲线图。
可扫任意字段（`faultline`、`Pm`、`Tunit`…），`values` 每个元素是该字段一整组取值。

---

## 5. 可追溯性

每次 `run_bcu` 在 `results/` 存 `snapshot_<时间戳>.{mat,txt}`，txt 是可读的参数清单，
可直接引用进论文的复现实验附录。扫参结果存 `sweep_<时间戳>.{csv,mat}`。

---

## 6. 验证与自检脚本（`verify/`）

| 脚本 | 作用 |
|------|------|
| `verify/v0_baseline.m` | 端到端基准 + SEP/CUEP 残差 + 能量法保守性自检，存 `baseline_reduced.mat` |
| `verify/v1_unit.m` | 5 项单元测试（SEP 处右端=0、势能守恒基准、路径积分自洽、RK4 步长收敛…） |
| `verify/v2_energy.m` | 有损/无损能量守恒对照 + 积分器阶数判据 |

改动核心代码后，重跑这三个脚本即可回归验证是否引入偏差。

---

## 7. 常见问题

- **`Unrecognized function 'setup_bcu_paths'`**：`run()` 会临时切到脚本目录；`setup_bcu_paths`
  已把项目根加入 path 修复此问题。若仍出现，先 `cd` 到 `matlab_platform` 再运行。
- **OSQP 报错**：本机未编译 `osqp_mex`。`setup_bcu_paths` 已预热 `have_feature('osqp')`
  缓存为 0，规避 MATPOWER `genpath` 引入未编译接口导致的间歇性崩溃。潮流走 Newton 法不需 osqp。
- **中文路径下命令行输出乱码**：仅命令行管道的 GBK/UTF-8 显示问题，数值无影响；MATLAB GUI 内正常。
- **`clear` 陷阱**：`Cal_MM_CCT` 会 `clear` 工作区。这是配置用"函数"（`bcu_config`）而非工作区
  变量注入、扫参用 `appdata(0)` 保存状态的根本原因。自定义脚本时请注意此点。
