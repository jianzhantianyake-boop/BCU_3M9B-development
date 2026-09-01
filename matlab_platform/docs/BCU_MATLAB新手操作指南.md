# BCU_3M9B MATLAB 平台新手操作指南

## 1. 这份指南解决什么问题

本指南面向第一次接触本项目的使用者，目标是让你以**安全、可复现的单次实验**方式运行 BCU_3M9B 平台。它只说明原有模型的运行顺序、参数位置、图窗含义和排错边界；不会替换原模型，也不会自动扫参、自动保存图片或批量生成结论。

当前默认且最适合新手复现的是 **3 机 9 母线（`case9_v2`）**。项目中虽然保留了 `case39_modified` 和对应注释参数，但把模型切到 39 母线前，必须同时检查发电机数、`preset.m`、`preset.d`、`preset.Pmpu`、`preset.xd1`、`preset.Epu` 的长度和顺序；这不是“只改一行 Case”即可安全完成的操作。

> 重要边界：本环境曾出现 MATLAB R2024a 启动期的 `errors_warnings` 插件加载问题。因此，本文档给出的是代码层面的可运行顺序，不把静态检查或 Python 冒烟结果当作 MATLAB 原生端到端验证。请在你自己的 MATLAB 中实际运行并保存命令窗口记录。

## 2. 平台地图与实验选择

| 目标 | 推荐控制器模式 | 原始运行链 | 主要结果 |
|---|---|---|---|
| 网络约简模型的 CCT/CUEP | `"reduced_cct"` | `Cal_MM_CCT` | `postfault.CUEP_delta`、逃逸点/MGP 及相平面图 |
| 网络约简模型的数值轨迹与能量图 | `"reduced_numerical"` | `Cal_MM_CCT` → `NumSim_MM_Gridframe` | `IterData`、转子角/角速度/能量图 |
| 网络约简模型二维稳定区域 | `"reduced_region"` | `Cal_MM_Static` → `Statable_Region` | `ep_set`、稳定区域图 |
| 结构保持模型的 CCT/CUEP | `"spm_cct"` | `Cal_MM_CCT_SPM` | 结构保持模型的 CUEP、MGP 与图窗 |
| 结构保持模型数值轨迹 | `"spm_numerical"` | `Cal_MM_CCT_SPM` → `NumSim_MM_Gridframe_SPM` | DAE 轨迹、节点角/电压/能量图 |
| 结构保持模型稳定区域 | `"spm_region"` | `Cal_MM_Static_SPM` → `Statable_Region_SPM` | 扩展平衡点集合与图窗 |
| 独立两机三维示例 | `"two_machine_region_3d"` | `Statable_Region_3D` | 三维平衡点/稳定区域探索图 |
| 独立 GFL 相关三维示例 | `"two_machine_region_3d_gfl"` | `Statable_Region_3D_GFL` | 原脚本定义的 GFL 相关探索图 |

`Statable_Region_3D_GFL.m` 中的 “GFL” 是原始文件命名。它应被理解为特定参数下的探索性示例，而不是对所有 GFL 工况均已完成验证的结论。

## 3. 第一次运行前的检查

1. 在 MATLAB 中将当前文件夹切换到项目根目录：`BCU_3M9B-main`。
2. 在命令窗口运行：

   ```matlab
   clear; close all;
   paths = setup_bcu_paths();
   which fsolve
   which ode78
   ```

3. `which fsolve` 应指向 MATLAB Optimization Toolbox；`which ode78` 应能找到项目或 MATLAB 可用实现。任何一个为空时，先解决环境问题，不要继续解释模型结果。
4. 如 MATLAB 本身启动失败或报 settings/plugin 错误，先记录完整报错、MATLAB 版本和 `ver` 输出。不要改动模型代码来绕过 MATLAB 启动问题。

## 4. 最推荐的运行方式：单次实验控制脚本

打开项目根目录的 [`run_bcu_beginner.m`](../run_bcu_beginner.m)，只修改：

```matlab
EXPERIMENT_MODE = "reduced_cct";
```

然后点击 MATLAB 编辑器的“运行”。控制脚本会：

1. 自动定位项目根目录；
2. 调用 `setup_bcu_paths` 配置 `B3_MM` 和项目自有 MATPOWER 接口路径；
3. 仅运行一个模式对应的原始脚本链；
4. 保留原始脚本生成的 base workspace 变量和图窗；
5. 不改模型参数、不扫参、不自动导出图件。

每次改模式前，请先执行：

```matlab
clear; close all;
```

这是必要步骤，不是形式步骤。原始项目有许多脚本通过 MATLAB **base workspace** 传递 `preset`、`prefault`、`fault`、`postfault` 和 CUEP 结果；混用上一次实验残留变量会使结果不可解释。

### 不要做的事

- 不要把 `run_bcu_beginner.m` 改成函数。它需要与原始脚本共享 base workspace。
- 不要在 CCT 与数值轨迹之间插入 `clear`、`clearvars` 或手动改写 `postfault`。
- 不要直接运行 `NumSim_MM_Gridframe.m` 或 `NumSim_MM_Gridframe_SPM.m` 作为首次实验。它们依赖 CCT 计算生成的 `postfault.CUEP_delta` 等数据。
- 不要同时打开两个控制脚本窗口或在两个 MATLAB 会话写同一输出变量后比较图窗。

## 5. 不使用控制器时的手动运行顺序

如需逐文件学习，可按下面顺序在 MATLAB 命令窗口执行。所有命令都假设当前文件夹是项目根目录。

### 5.1 网络约简模型

**只计算 CCT / CUEP：**

```matlab
clear; close all;
setup_bcu_paths();
run(fullfile('B3_MM', 'Cal_MM_CCT.m'));
```

**计算 CCT 后再做数值轨迹：**

```matlab
clear; close all;
setup_bcu_paths();
run(fullfile('B3_MM', 'Cal_MM_CCT.m'));
run(fullfile('B3_MM', 'NumSim_MM_Gridframe.m'));
```

**只做二维稳定区域搜索：**

```matlab
clear; close all;
setup_bcu_paths();
run(fullfile('B3_MM', 'Cal_MM_Static.m'));
run(fullfile('B3_MM', 'Statable_Region.m'));
```

### 5.2 结构保持模型（SPM）

**只计算 CCT / CUEP：**

```matlab
clear; close all;
setup_bcu_paths();
run(fullfile('B3_MM', 'Cal_MM_CCT_SPM.m'));
```

**计算 CCT 后再做数值轨迹：**

```matlab
clear; close all;
setup_bcu_paths();
run(fullfile('B3_MM', 'Cal_MM_CCT_SPM.m'));
run(fullfile('B3_MM', 'NumSim_MM_Gridframe_SPM.m'));
```

**稳定区域搜索：**

```matlab
clear; close all;
setup_bcu_paths();
run(fullfile('B3_MM', 'Cal_MM_Static_SPM.m'));
run(fullfile('B3_MM', 'Statable_Region_SPM.m'));
```

### 5.3 独立两机三维示例

这两份脚本有各自的 `Z1`、`Z2`、`Zl`、`Pm1`、`Pm2`、`H1`、`H2`、`D1` 和 `D2`；它们不使用 9 母线初始化链。

```matlab
clear; close all;
setup_bcu_paths();
run(fullfile('B3_MM', 'Statable_Region_3D.m'));
% 或运行：Statable_Region_3D_GFL.m
```

## 6. 参数应改哪里，以及哪些地方暂时不要碰

修改前先复制原始数值到实验记录中。一次只改一类参数，运行后确认图窗和命令窗口正常，再进行下一次修改。

| 实验意图 | 修改文件与变量 | 单位 / 约束 | 新手提醒 |
|---|---|---|---|
| 切换约简/结构保持模型的计算策略 | `Cal_MM_Static.m` 或 `Cal_MM_Static_SPM.m` 中 `preset.PathEnergyCal`、`preset.EquCal` | `PathEnergyCal` 为离散策略标志；无物理单位 | 这是算法设置，不等于物理故障参数；与数值脚本内的能量对比设置区分开。 |
| 发电机惯性、阻尼 | 同一初始化文件中的 `preset.m`、`preset.d` | 原始模型的标幺/归一化参数；请同时记录 `Basevalue.omegab` | 3 个数值必须与 9 母线 3 台发电机的顺序一致。 |
| 机械功率、暂态电抗、内部电势 | `preset.Pmpu`、`preset.xd1`、`preset.Epu` | pu | 四个发电机向量的长度必须一致；不要只改其中一个向量的长度。 |
| SPM ZIP 负荷比例 | `Cal_MM_Static_SPM.m` 的 `preset.PloadZIP`、`preset.QloadZIP` | 比例，无单位；顺序为 Z/I/P | 改完后应检查比例定义及是否仍满足你的研究设定。 |
| 选择默认 9 母线工况 | 初始化文件的 `Case=case9_v2;` | 无 | 这是当前新手路线。切到 `case39_modified` 前需重配全部发电机参数和维度。 |
| 故障线路与故障端 | 初始化文件的 `preset.faultline=[9;6]`、`preset.faultposition` | 母线编号；`faultposition` 为线路端索引 | 保持 `faultline` 为两端母线的列向量。先确认线路确实存在于所选 case。 |
| CCT 脚本内的故障持续时长与步长 | `Cal_MM_CCT.m` 的 `Tfault`、`Tunit`；或 SPM 对应文件 | s | `Tfault` 是用于该 CCT 过程的设定；减小 `Tunit` 会显著增加计算量。 |
| 数值轨迹的故障时刻、切除时刻、总时长、步长 | `NumSim_MM_Gridframe*.m` 的 `Iter.Tfault`、`Iter.Trecover`、`Iter.Ttotal`、`Iter.Tunit` | s | 必须满足 `0 ≤ Tfault < Trecover ≤ Ttotal`。故障持续时间为 `Trecover - Tfault`。 |
| 数值轨迹的能量项比较 | 两个 `NumSim_MM_Gridframe*.m` 中的 `preset.PathEnergyCal` 赋值段 | 算法标志，无物理单位 | 原脚本会分别使用 0、1、10 等策略作比较；先理解对应图例，再删改比较支路。 |

### 参数修改的最小安全流程

1. 在 `docs/实验记录模板.md`（若你自行创建）或笔记中记录原始值和新值。
2. 只修改一个文件中的一组同类变量，例如只修改 `Iter.Trecover`。
3. 运行对应的 `*_cct` 或 `*_numerical` 单次模式。
4. 检查是否出现维度错误、`fsolve` 未收敛、负电压或 CUEP 缺失。
5. 保存图窗、工作区变量名和命令窗口输出；再开始下一组设置。

## 7. 图窗怎么得到、怎么看、怎么手动保存

### 7.1 CCT / CUEP 图窗

运行 `reduced_cct` 或 `spm_cct` 后，原脚本会绘制以第 2、3 台机相角为坐标的相平面信息。通常包括：

- 故障轨迹；
- 逃逸点（escape point）；
- 故障前 SEP；
- MGP；
- `postfault.CUEP_delta` 对应的 CUEP 标记。

这些图用于核对原始脚本中各状态点的相对位置；不能仅凭某个标记点存在就推断所有扰动下都稳定。

### 7.2 网络约简数值轨迹图

`reduced_numerical` 会绘制随时间变化的：

- 转子/相对相角；
- 角速度与 COI 相关量；
- 多种势能、动能、阻尼能量及其差异；
- 相平面轨迹。

横轴通常是时间 s；相角通常为 rad；速度通常为 rad/s 或相对归一化量。实际变量名以 `IterData` 的字段为准。建议先在命令窗口运行 `whos IterData`，再读取需要的字段。

### 7.3 SPM 数值轨迹图

`spm_numerical` 除发电机角度和角速度外，还会绘制网络节点角度、网络节点电压，以及扩展能量项。SPM 的状态量更多，图窗数也更多；先确认每张图的图例和脚本附近的 `ylabel`，不要把网络节点电压曲线误认为发电机内部电势。

### 7.4 独立绘图脚本

- `Plot_3Dstate.m` 不是主入口。它会 `load('Data_InitState_Dp20_H0.1_ZL1.mat')`，仅当该 MAT 文件存在于当前工作目录且包含脚本所需变量时才可运行。先在 MATLAB 中执行 `isfile('Data_InitState_Dp20_H0.1_ZL1.mat')`；返回 `0` 时不要运行它。
- `plottmp.m` 依赖 base workspace 中已有的 `Group.delta_stb` 和 `Group.CUEP`。只有你已生成并检查过 `Group` 结构体时才可运行；它不是数据生成脚本。
- `vectorfield_cal.m` 是函数，不是按钮式入口。其调用形式为：

  ```matlab
  f = vectorfield_cal(thetac, Yred_post, preset);
  ```

  其中 `thetac` 为角度状态（rad），`Yred_post` 为故障后约简导纳矩阵（pu 相关），`preset` 为初始化产生的参数结构体。

### 7.5 手动保存一张已确认的图

本平台的控制脚本不会自动导出。检查图窗正确后，可在命令窗口手动执行：

```matlab
exportgraphics(gcf, '我的实验图.png', 'Resolution', 300);
```

建议文件名包含模型、故障线路、清除时间和日期，例如 `reduced_9-6_tclear0p24_2026-08-27.png`。保存前请先确认当前活动图窗就是你要保存的那一张。

## 8. 常见错误与处理顺序

| 现象 | 常见原因 | 首先做什么 |
|---|---|---|
| `Unrecognized field name "CUEP_delta"` | 直接跑数值轨迹，或前一轮 CCT 未找到 CUEP | `clear; close all;` 后先完整运行同一模型的 CCT 模式；确认命令窗口出现 CUEP 相关结果，再跑轨迹。不要用全零向量伪造 CUEP。 |
| `Undefined function 'fsolve'` | Optimization Toolbox 不可用或路径异常 | 运行 `ver` 和 `which fsolve`，恢复工具箱后再运行。 |
| 找不到 `ode78` | 自定义求解器未在路径中或文件缺失 | 先运行 `setup_bcu_paths()`，再执行 `which ode78`；不要自行替换为其他求解器并声称结果等价。 |
| 发电机数不匹配 | case 文件与 `preset` 向量长度/顺序不一致 | 回到初始化脚本，核查 `preset.m`、`preset.d`、`preset.Pmpu`、`preset.xd1`、`preset.Epu` 和 case 的发电机数量。 |
| 稳定区域脚本极慢 | 网格点多、`fsolve` 最大迭代次数高 | 保留一份原始参数，先在独立副本中降低搜索网格密度用于学习；不要把低分辨率图当作正式结论。 |
| 图窗为空或数据不连续 | 初始化没完成、工作区被清空，或求解器未收敛 | 查看命令窗口的第一个报错；重新按完整链执行。不要只重跑最后一个绘图脚本。 |
| MATLAB 启动/设置插件报错 | MATLAB 本体环境问题 | 记录报错、MATLAB 版本和 `ver`；先修复 MATLAB，再谈模型验证。 |

## 9. 每次实验应记录的最小信息

至少记录：

1. MATLAB 版本、操作系统、`ver` 输出和运行日期；
2. 实验模式与入口文件；
3. case 文件、发电机参数、故障线路、故障/切除时刻、积分步长；
4. 是否找到 CUEP，以及对应的命令窗口提示；
5. 生成的图文件名、关键工作区变量和任何警告/报错；
6. 结果是“实际 MATLAB 运行通过”“近似接口结果”还是“仅静态检查”。

这样才能把“看到一张图”与“得到可复现研究结论”区分开来。

## 10. 本次新增控制层的边界

`run_bcu_beginner.m` 只是便于新手进入原有实验平台的单次调度器。它不改变：

- 原始 BCU/能量函数/稳定判据；
- 原始模型参数与默认 case；
- 原始绘图数据和图形风格；
- MATPOWER 官方源代码；
- MATLAB 与 Python 结果之间尚未建立的严格定量交叉验证。

若你准备做参数扫描、批量导出、39 母线重构、控制器对比或论文级统计，请先另行定义实验设计、验收标准和结果保存规范，再在独立任务中实施。
