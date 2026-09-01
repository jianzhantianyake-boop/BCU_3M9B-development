# BCU_3M9B 学习与运行指南

## 1. 这套代码要解决什么问题

本项目研究的是三机九母线系统受到故障后，如何判断故障切除得是否足够及时。
代码同时提供两类思路：

1. 直接数值积分：分别模拟预故障、故障和故障后网络，观察转子角、转速和能量；
2. BCU / Energy Function：寻找故障轨迹退出点，追踪稳定边界，选取 MGP，求解
   CUEP，再用临界能量估计临界切除时间 CCT。

这里的“稳定”必须绑定到具体模型、故障、时间窗和判据。能量法得到的 CCT 是
能量判据结果，时域仿真的 CCT 是数值轨迹判据结果，二者需要交叉比较，不能在
没有验证时互相冒充。

## 2. 运行前准备

推荐环境：

- MATLAB R2024a 或更高版本；
- 本仓库随附的 MATPOWER 7.1；
- `fsolve`、`optimoptions` 等函数需要 MATLAB Optimization Toolbox；
- 部分入口使用 `ode78`，运行前应检查当前 MATLAB 版本是否提供该求解器；
- Windows 路径可以包含中文，但建议不要移动 `B3_MM` 与 `C1_Matpower` 的相对层级。

项目不再依赖原作者的 OneDrive 目录。路径由项目根目录下的
`setup_bcu_paths.m` 自动配置。

## 3. 第一次启动

在 MATLAB 当前目录切换到 `BCU_3M9B-main` 根目录后执行：

```matlab
paths = setup_bcu_paths();
which runpf -all
which case9_v2 -all
which Fun_ResultBack -all
```

三个 `which` 命令都应指向当前仓库内的文件。若 MATLAB 当前目录不是项目根目录，
也可以传入绝对路径：

```matlab
paths = setup_bcu_paths('C:\path\to\BCU_3M9B-main');
```

旧入口 `Y.m` 仍可直接运行：

```matlab
Y
```

但它只是 `setup_bcu_paths()` 的兼容包装，不再负责切换当前目录或添加作者电脑
上的路径。

## 4. 推荐学习顺序

### 第一步：两机模型

先阅读：

- `B3_MM/f_2m.m`：两机完整状态；
- `B3_MM/f_2m_fault.m`：故障期间状态方程；
- `B3_MM/f_2m_reduce.m`：消去共同速度后的二阶模型。

学习目标：理解 `delta12`、`omega12`、`omegasum`、机械功率 `Pm`、
电磁功率 `Pe`、惯性和阻尼项之间的关系。两机模型是理解多机 COI 坐标和
能量函数的最小入口。

### 第二步：9-bus 系统和潮流桥接

推荐阅读：

1. `C1_Matpower/matpower7.1/data/case9_v2.m`；
2. `C1_Matpower/matpower7.1/Fun_ResultBack.m`；
3. `B3_MM/Fun_Cal_GenEMF.m`；
4. `B3_MM/Fun_RXB2Yfull.m`。

当前默认系统有 3 台发电机，位于 1、2、3 号母线；默认故障线路为 `[9; 6]`，
故障位置参数为 0，因此故障母线为 9 号母线。切换 case 或故障位置时，必须
同步检查发电机数量、`preset.Pmpu`、`preset.m`、`preset.d`、`preset.Epu`
和矩阵维度。

### 第三步：Y-bus 与网络约简

阅读：

- `B3_MM/Fun_Yfull2Yfull.m`：结构保持模型的节点重排；
- `B3_MM/Fun_Yfull2Yred.m`：Kron 消元；
- `B3_MM/Cal_MM_Static.m`：网络约简初始化；
- `B3_MM/Cal_MM_Static_SPM.m`：结构保持初始化。

约简模型把网络节点消去，得到发电机之间的 `Yred`；结构保持模型保留网络节点
角度和电压，后续使用更大的 DAE 状态。两者不是同一组状态变量，不能只比较
数组长度或单个角度就宣称模型等价。

### 第四步：SEP、CUEP 和 MGP

推荐顺序：

1. `Fun_SEPfslove.m` / `Fun_SEPfslove_SPM.m`；
2. `Fun_SEPiteration.m` / `Fun_SEPiteration_SPM.m`；
3. `Fun_SEPcheck.m`；
4. `Fun_Cal_MGP.m` / `Fun_Cal_MGP_SPM.m`；
5. `Fun_Cal_MGP_singletraj.m` / `Fun_Cal_MGP_singletraj_SPM.m`。

SEP 是平衡点。MGP 是沿稳定边界追踪时用于初始化 CUEP 搜索的候选点。
CUEP 是临界不稳定平衡点候选。`Fun_SEPcheck.m` 的残差检查是必要条件，
但残差小本身不能证明点一定是正确的 CUEP 或全局临界点。

### 第五步：能量函数

阅读：

- `Fun_Cal_PotentialEnergy.m` / `_SPM.m`；
- `Fun_Cal_DampingEnergy.m`；
- `Fun_Cal_CCT_Energy.m` / `_SPM.m`；
- `Fun_AEiteration_SPM.m`。

重点追踪：

- 动能 `Ek` 如何由转速偏差得到；
- 势能 `Ep` 如何沿功角/网络状态积分；
- 阻尼能如何累计；
- CUEP 处的临界能量如何与故障轨迹能量比较。

能量的数值大小依赖参考平衡点、坐标、单位和积分近似方式。Ray approximation
与分段梯形积分的结果应分别记录。

### 第六步：CCT 与时域仿真

网络约简模型：

```matlab
paths = setup_bcu_paths();
Cal_MM_CCT
```

该入口会依次执行静态初始化、故障轨迹退出点、MGP、CUEP 和能量 CCT，并绘制
角度平面上的关键点。成功时应在命令行看到 Exit point、MGP 和 CUEP 等信息。

结构保持模型：

```matlab
Cal_MM_CCT_SPM
```

数值时域仿真应先完成对应静态初始化。网络约简模型：

```matlab
Cal_MM_Static
NumSim_MM_Gridframe
```

结构保持模型：

```matlab
Cal_MM_Static_SPM
NumSim_MM_Gridframe_SPM
```

仿真由预故障、故障、故障后三段组成。先检查 `Iter.Tfault`、`Iter.Trecover`、
`Iter.Ttotal` 和 `Iter.Tunit`，再解释图形。离散 RK4、`ode78` 和结构保持
`ode15s` 的结果应保留各自的求解器和容差信息。

## 5. 核心变量和单位

| 变量 | 含义 | 常见单位/约定 |
|---|---|---|
| `delta` / `theta` | 发电机或网络节点相角 | rad |
| `omega` | 发电机角速度 | rad/s |
| `omegac` | 相对于 COI 的速度 | rad/s |
| `Epu`, `V_net` | 发电机内部电势、网络电压幅值 | pu |
| `Pm`, `Pe` | 机械输入功率、电磁输出功率 | pu |
| `Yfull`, `Yred` | 完整/约简复导纳矩阵 | pu 导纳 |
| `SEP_delta`, `CUEP_delta` | 平衡点相角 | rad |
| `Iter.T*`, `Tunit` | 仿真时间和离散步长 | s |
| `m`, `d` | 原始代码中的惯量/阻尼系数 | 按模型标幺/归一化定义 |

不要仅根据变量名把 `m` 直接当作 SI 制惯量常数；应结合
`Basevalue.omegab` 和对应摆动方程判断其归一化方式。

## 6. MATLAB 脚本依赖的注意事项

许多原始函数通过 `evalin('base', ...)` 读取 `preset`、`prefault`、`fault`、
`postfault`、`Basevalue` 和 `system`。这意味着：

- 直接单独调用某个函数通常会报“变量不存在”；
- 应先运行对应的 `Cal_MM_Static*.m` 或在测试中显式构造 base workspace；
- 不要在静态初始化和后续仿真之间随意执行 `clear`；
- `NumSim_MM_Gridframe*.m` 中部分初始化代码被注释掉，必须按本指南先运行静态脚本；
- 图形窗口和命令行输出是当前原始实现的主要结果载体，暂不等同于结构化结果文件。

工具箱/函数排查命令：

```matlab
which fsolve -all
which optimoptions -all
which ode78 -all
ver
```

## 7. 当前验证边界

本阶段已经完成项目路径去硬编码、教学注释和文档化。MATLAB R2024a 已安装，
但当前机器的命令行启动仍报告：

```text
failed to load settings errors_warnings plugin
```

在该 MATLAB 环境问题解决前，不能把 `Cal_MM_Static`、CCT 或暂态仿真标记为
原生 E2E 通过。后续验证应保存 MATLAB 版本、工具箱、命令、命令行输出和结果
图形，并明确区分静态检查、脚本启动、数值运行和结果交叉验证。

## 8. 后续学习与研究工作

本阶段没有实现 Python 版本，也没有宣称 MATLAB/Python 数值一致。后续应先定义
稳定的输入/输出数据结构，再分别实现约简模型和结构保持模型，并使用同一工况、
故障、时间网格、容差和判据进行定量交叉验证。
