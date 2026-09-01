# 工作日志：BCU_3M9B 可靠性验证与操作入口平台化

- 日期：2026-08-31
- 环境：Windows 11 + MATLAB R2024a Update 8 + Optimization Toolbox + VSCode（`MathWorks.language-matlab`）
- 代码来源：本平台为 [yifanz1125/BCU_3M9B](https://github.com/yifanz1125/BCU_3M9B)（Imperial College，作者 yz7521）的可移植性改写版
- 目标：(1) 验证代码可靠性并定位潜在错误；(2) 搭建配置驱动的操作入口与批量扫参；(3) 版本控制

---

## 1. 环境配置

- 确认 MATLAB R2024a 在 `C:\Program Files\MATLAB\R2024a`，`matlab -batch` 可命令行驱动
- 确认 Optimization Toolbox / `fsolve` / `ode78` 可用
- 写入 `.vscode/settings.json`（`matlab.installPath`、连接时机、诊断开关）
- VSCode 内 ▶ Run / `Ctrl+Enter` 逐节运行 / 断点调试均可用

---

## 2. 可靠性验证方案（四层）与结果

采用"从便宜到昂贵"的分层验证。全部脚本在 `verify/`，可复现。

### 第 0 层：静态审查 + 冒烟
- 端到端跑通 `Cal_MM_Static + Cal_MM_CCT`（约 55 s）
- 识别死代码：`Fun_Cal_PotentialEnergy` 中 `i=preset.m; d=preset.d;` 无用赋值、`Ep3_rad` 段算而不用

### 第 1 层：单元级解析解/守恒律（`verify/v1_unit.m`）——**5/5 全 PASS**
| 测试 | 结果 |
|------|------|
| U1 摆动方程右端在 SEP 处 = 0 | 1.6e-16 |
| U2 势能起点=终点 → 三项 = 0 | 0 |
| U3 无损网络 → 路径项 Ep3 = 0 | 0 |
| U4 路径积分自洽（Ray vs 200 段梯形） | 相对差 0.00% |
| U5 RK4 步长减半，退出点收敛 | 差 5.6e-5 rad |

### 第 2 层：内部一致性 / 守恒律（`verify/v0_baseline.m`, `v2_energy.m`）
- SEP/CUEP 平衡残差：约 1e-16（PASS）
- CUEP 显著区别于 SEP：|CUEP−SEP| = 2.4543（PASS）
- 能量法保守性：LEA-CCT(0.2275) ≤ REA-CCT(0.2434)，差 6.5%（PASS）
- 能量守恒：无损网络下总能量漂移随步长 **严格一阶（O(h)，比值 2.00/2.00）** 下降 → 能量函数实现正确

### 第 4 层：与上游原版对照
- 克隆 `yifanz1125/BCU_3M9B`，逐文件对比 `B3_MM`（去注释后）
- **所有核心数值函数零实质差异**（SEP/CUEP 求解、能量函数、CCT 判据、积分器一行未改）
- 仅 4 个文件有实质改动，均为可移植性/兼容性：
  - 新增 `setup_bcu_paths()`；`roundn(x,-3)` → `round(x,3)`（等价）
  - 硬编码 `C:\Users\yz7521\OneDrive...` → `projectPaths.*`（去绝对路径）
  - `Statable_Region_3D` 绘图代码移位 + 补 `maxabs` 局部函数

**结论：改写是纯可移植性改造，数值正确性完整继承自上游，且无回归。**

---

## 3. 发现的问题（均非致命）

1. **积分器精度缺陷（唯一实质缺陷，已量化）**：`Fun_Cal_Exitpoint` 名为 RK4，但四个 stage
   未在中间点重算电磁功率 `Pe`；`d=0` 时 `ω` 分量退化为前向欧拉一阶。经 `v2_energy.m`
   步长阶数判据证实（漂移比值精确 = 2.00）。
   - 影响：长时间保守积分有 O(h) 能量漂移；实用 CCT 场景（Tunit=1e-4、故障窗 ~0.24 s）
     影响可忽略（U5：步长减半退出点仅变 5.6e-5 rad）。
   - 建议：若做长轨迹能量分析，改为真正耦合 RK4 或辛积分器。
2. **死代码**：`Fun_Cal_PotentialEnergy` 的 `i=preset.m; d=preset.d;`、`Ep3_rad` 段。
3. **方法固有近似**：能量法 CCT 比真实 CCT 保守 6.5%，为 TEF 有损网络路径近似的正常表现。

---

## 4. 操作入口平台化改造

原入口 `run_bcu_beginner.m` 只能改一行 mode，参数散落 2–3 个文件且深藏 `Cal_MM_Static`。
本次做配置驱动重构，**不改任何核心方程**。

### 新增文件
| 文件 | 作用 |
|------|------|
| `bcu_config.m` | 唯一参数配置入口（6 组参数，详注单位/范围/约束） |
| `bcu_pick.m` | 参数覆盖钩子：配置有值则覆盖，否则用原默认（保证零回归） |
| `bcu_override.m` | 运行时覆盖通道（persistent，供扫参用，`clear` 不清） |
| `bcu_validate_config.m` | 运行前校验（维度/支路存在性/步长…），fail-fast |
| `run_bcu.m` | 主入口：加载→校验→摘要→快照→跑链路→自检 |
| `run_bcu_sweep.m` | 批量扫参驱动，输出 CSV + CCT 曲线 |
| `docs/BCU平台操作入口指南.md` | 使用文档 |

### 改动文件（最小侵入）
- `Cal_MM_Static.m`：8 处硬编码参数包 `bcu_pick(BCUCFG,...)`，默认值原样保留
- `Cal_MM_CCT.m`：`Tfault/Tunit` 走配置
- `setup_bcu_paths.m`：补 `addpath(项目根)`（修 `run()` 切目录后根脚本失联）+ 预热
  `have_feature('osqp')`（修 MATPOWER `genpath` 引入未编译 osqp 接口导致的间歇崩溃）

### 关键技术难点与破解
- **`Cal_MM_CCT` 的 `clear` 会清空工作区**：配置用"函数"而非工作区变量注入；扫参用
  `appdata(0)` 保存循环状态、`persistent` 保存参数覆盖，两者都不受 `clear` 影响，
  循环得以跨多次 `clear` 稳定推进。
- **`run()` 临时切目录**：靠 `setup_bcu_paths` 把项目根加入 path 解决。
- **OSQP 间歇崩溃**：`setup_bcu_paths` 末尾预热探测缓存。

### 无回归 + 配置生效验证
- 默认配置下 LEA/REA-CCT、CUEP、残差与改造前 **bit 级一致**（零回归）
- 阻尼比 0.1→0.3 时 CCT 0.2275→0.2308（↑，物理正确），证明配置真正驱动结果
- 扫参 4 组（damping 0.05/0.1/0.2/0.3）CCT 严格单调，CSV/图正常输出

---

## 5. 文件与目录

```
matlab_platform/
├── bcu_config.m            ← 你要编辑的参数文件
├── bcu_pick.m / bcu_override.m / bcu_validate_config.m
├── run_bcu.m               ← 单次实验主入口
├── run_bcu_sweep.m         ← 批量扫参
├── run_bcu_beginner.m      ← 旧入口（保留，兼容）
├── setup_bcu_paths.m
├── B3_MM/                  ← 核心算法（仅 Static/CCT 加了配置注入）
├── verify/                 ← v0/v1/v2 验证脚本 + baseline_reduced.mat
├── results/                ← 运行快照与扫参输出（git 忽略）
├── docs/                   ← 操作指南等
└── C1_Matpower/            ← MATPOWER 7.1（外部依赖，git 忽略）
```

---

## 6. 下一步方案

1. **SPM 链路参数注入**：将 `bcu_config` 覆盖接入 `Cal_MM_Static_SPM.m`（结构保持模型），
   使 SPM 链路与 reduced 链路共享同一套参数。
2. **修积分器**：把 `Fun_Cal_Exitpoint` 改为真正的耦合 RK4（中间 stage 重算 `Pe`），
   用 `v2_energy.m` 验证漂移降到 O(h⁴)。
3. **SPM 链路验证**：为 SPM 补一套 v0/v1/v2 同类验证（目前只验了约简模型）。
4. **教科书基准对照**：与 Anderson-Fouad / Sauer-Pai 的 WSCC 3 机 9 节点公开 CCT 对齐。
5. **清理死代码**：移除 `Ep3_rad` 段与无用赋值（先用 verify 回归确认无影响）。
6. **扫参增强**：支持二维扫参（如 faultline × damping 网格）与并行。
