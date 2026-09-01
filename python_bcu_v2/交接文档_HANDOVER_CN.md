# BCU 暂态稳定平台 · 开发交接文档

> 目的：让接手者在 30 分钟内看懂整个项目、跑通、并能安全地继续开发。
> 最后更新：2026-08-31。单位约定：角度 rad，角速度 rad/s，时间 s，功率/电压/导纳 pu，频率 Hz。

---

## 0. 30 秒速览

- **做什么**：电力系统**暂态稳定（功角首摆）+ 直接法（BCU / 能量函数）**的仿真研究平台。
- **三个组成**（互不干扰）：
  1. `matlab_platform/` —— **已验证的 MATLAB 平台**（参考真值来源，**只读，禁止改**）。
  2. `python_bcu/` —— **v1**：MATLAB→Python 的忠实复写，**已冻结**为可用基线。
  3. `python_bcu_v2/`（本目录）—— **v2**：在 v1 之上做正确性验证与数值/建模增强，**只增量、不改 v1**。
- **可信度已建立**：v2 与 MATLAB 逐层交叉验证（T3）核心层吻合到 **1e-10~1e-11**；并有不依赖任
  何实现的物理不变量 + 等面积闭式解金标准（P0）。
- **两个私有 GitHub 库**：`python_platform_lab`（v1）、`python_platform_lab_v2`（v2）。

## 1. 最快上手（接手者照这个顺序做）

```bash
# 1) 进 v2 目录
cd python_bcu_v2

# 2) 跑一遍全部验证, 确认环境与代码健康
python run_validation.py     # P0: 物理不变量 + 等面积金标准, 应 8/8
python run_matlab_xval.py    # T3: 与 MATLAB 逐层对比, 应 8/8
python test_config.py        # 配置系统 + T3 回归, 应 5/5

# 3) 用命令行入口跑实验(这是主操作界面)
python cli.py show           # 看当前配置
python cli.py list           # 看可选 mode / case
python cli.py run            # 用 config.yaml 跑(含运行后自检)

# 4) 调参: 改 config.yaml, 或命令行临时覆盖
python cli.py run --case case39_modified --auto-params   # 一键切 39 母线
python cli.py run --mode reduced_region --grid 15
```

若上面三个验证都通过，说明你已经把整套东西跑起来了。**下一步开发看第 7、8 节。**

## 2. 环境与依赖

| 项 | 值 / 位置 |
|---|---|
| OS | Windows 11 |
| Python | 3.12.1（`C:\Users\WangS\AppData\Local\Programs\Python\Python312`） |
| 必需库 | `numpy` 2.2.4、`scipy` 1.17.1 |
| 可选库 | `matplotlib` 3.10.1（画图）；**`PyYAML` 未装**（config.py 内置极简解析已兜底） |
| MATLAB | R2024a（`C:\Program Files\MATLAB\R2024a`，`matlab` 不在 PATH；T3 用已导出的 .mat，**无需运行 MATLAB**） |
| Git/gh | gh 2.98 已装并登录（账号 jianzhantianyake-boop，keyring） |

**注意**：v2 通过 `bcu_v2/__init__.py` 把兄弟目录 `../python_bcu` 加入 `sys.path`，因此能
`import bcu_3m9b`。跑 v2 脚本请在 `python_bcu_v2/` 目录下。

## 3. 三个平台的关系与目录

```
BCU_3M9B-main/
├── matlab_platform/        # 【只读参考】已验证 MATLAB 平台
│   ├── run_bcu.m           #   主入口(配置驱动)
│   ├── bcu_config.m        #   唯一配置文件(v2 的 config.yaml/cli 就是镜像它)
│   ├── B3_MM/              #   核心 .m(Cal_MM_Static/CCT, Fun_*, Statable_Region* ...)
│   └── verify/baseline_reduced.mat   # ★ T3 参考真值(Yred/SEP/CUEP/CCT...)
├── python_bcu/             # 【v1, 冻结】-> 私有库 python_platform_lab
│   └── bcu_3m9b/           #   核心包(见 4.1); main.py 8 模式控制台; scripts/ tests/
└── python_bcu_v2/          # 【v2, 本目录】-> 私有库 python_platform_lab_v2
    ├── bcu_v2/             #   增强包(见 4.2)
    ├── config.yaml         #   ★ 唯一需要编辑的参数文件
    ├── cli.py              #   ★ 命令行入口(run/show/list/validate/xval)
    ├── run_validation.py run_matlab_xval.py   # 一键验证
    ├── test_*.py           #   各阶段验证脚本
    └── *.md                #   README / 路线图 / 本交接文档
```

## 4. 模块地图（每个文件干什么）

### 4.1 v1 核心包 `python_bcu/bcu_3m9b/`（冻结，勿改）

| 模块 | 职责 |
|---|---|
| `types.py` | 全部数据结构（CaseData/PFData/Preset/BaseValue/NetworkState/Trajectory/StaticResult/CCTResult） |
| `cases.py` | 内置案例 `case9_v2()`、`case39_modified()` |
| `powerflow.py` | 交流潮流 `solve_power_flow`、`to_pfdata`、发电机内电势 |
| `network.py` | 导纳组装、故障设置、Kron 约简 `kron_reduce` |
| `equilibrium.py` | SEP/CUEP 残差与求解 `solve_sep`、`electrical_power` |
| `dynamics.py` | 暂态轨迹 `integrate_reduced`、退出点、时域 CCT |
| `energy.py` | 势能/动能/MGP/能量法 CCT `energy_cct`、`find_mgp` |
| `spm.py` | 结构保持模型（原始逐步牛顿版，偶发不收敛） |
| `two_machine.py` | 两机模型 + 平衡点分类 |
| `stability_region.py` | 稳定域/向量场工具 |
| `numerics.py` | 自写有限差分牛顿 `newton_solve`、RK4 |
| `bcu.py` | 一键入口 `build_static_result`、`run_bcu_experiment`、`default_preset` |
| `experiments.py` | **8 个实验模式** `mode_*` + `MODES` 字典（reduced/spm × cct/numerical/region、two_machine 3D/GFL） |
| `matlab_compat.py` | 与 MATLAB `Fun_*` 同名薄封装 |
| `plotting.py` | 基础绘图 |

**v1 顶层**：`main.py`（中文菜单控制台，8 模式 + 工具项，`python main.py`）、`draw_plots.py`、
`平台操作指南_CN.md`、`scripts/run_*.py`、`tests/smoke_test.py`。

### 4.2 v2 增强包 `python_bcu_v2/bcu_v2/`

| 模块 | 阶段 | 职责 |
|---|---|---|
| `smib.py` | P1.2 | SMIB + **等面积准则闭式 CCT**（验证金标准） |
| `cct.py` | P1.2 | **事件驱动 + 二分的精确 CCT**（solve_ivp 事件；SMIB & 3 机约简） |
| `fixes.py` | P0.2 | v1 五隐患的修正版包装（见第 6 节） |
| `invariants.py` | P0.1 | T1 物理不变量 + T2 独立参照（SMIB/scipy）检查集 |
| `spm_dae.py` | P1.3 | **严格 DAE 级 SPM**：约束流形降阶 ODE + 连续法（根治代数解脆弱） |
| `solvers.py` | P1.1 | SciPy 求解器统一层（root / v1 牛顿可开关，solve_sep_scipy，solve_ivp） |
| `systems.py` | P2.1 | **通用装配层** `build_preset` + **39 母线** `case39_dynamic` |
| `loads.py` | P2.3 | ZIP 负荷模型 |
| `models.py` | P2.2/2.4 | 单机级参考模型：one-axis 同步机、GFM 下垂逆变器 |
| `matlab_xval.py` | T3 | 与 MATLAB 平台逐层交叉验证（读 baseline_reduced.mat） |
| `cuep.py` | 方向1 | **通用能量法 CUEP + LEA CCT**（任意 ngen 的 closest-UEP 法，39 母线也能出 LEA） |
| `config.py` | — | 集中配置：加载 config.yaml / 校验 / 建 static / 摘要 / 快照 |

**v2 顶层**：`config.yaml`（编辑入口）、`cli.py`（命令行）、`run_validation.py`、
`run_matlab_xval.py`、`test_{fixes,spm_dae,solvers,p2,config}.py`。

## 5. 开发历程（时间线，理解“为什么长这样”）

1. **MATLAB→Python 复写（v1 起源）**：把 `B3_MM/` 的 MATLAB 工程改写成显式 Python 数据结构
   （dataclass 替代 base workspace）。
2. **注释规范化**：docstring 改为**中文「操作说明」风格**（使用方法/参数/返回/步骤），**程序
   输出改英文**，代码注释用**半角标点**（用户编辑器会破坏全角标点）。
3. **中央控制台**：`main.py` 做成 8 模式菜单（对应 MATLAB EXPERIMENT_MODE）；补齐 reduced/SPM
   六模式 + 两机 3D/GFL（发现 GFL 就是 `f_2m` 换参数，非新模型）。
4. **v2 · P0 正确性**：验证套件（物理不变量 + 等面积金标准）；修 v1 五隐患。
5. **v2 · P1.2**：事件驱动精确 CCT（发现 v1 网格 REA=0 是 bug）。
6. **v2 · P1.3**：严格 DAE 级 SPM（连续法根治代数解脆弱，DAE 7/7 vs v1 3/7）。
7. **v2 · P1.1 + P2**：scipy 求解器层；通用装配层（9→39 母线）；ZIP 负荷；单机 one-axis/GFM。
8. **v2 · T3 + 配置/CLI**：与 MATLAB 逐层交叉验证（8/8，并定位 v1 的 MGP 缺陷）；集中配置 +
   命令行入口（镜像 MATLAB bcu_config.m/run_bcu.m）。

## 6. 关键设计决策与约定（改代码前必读）

- ~~**v1 冻结**~~ → **v1 冻结已于 2026-09-01 由负责人解除**：6 个已知隐患已就地回灌修复(见第 7 节),
  v1 入口也移植成 v2 模式(`python_bcu/config.yaml` + `bcu_3m9b/config.py` + `cli.py`，保留 `main.py` 菜单)。
  **v1 现依赖 scipy**(原纯 numpy)。今后改 v1 仍需先跟负责人确认;改 v1 后**必须重跑 v2 全套**(v2 `import
  bcu_3m9b`, 会随 v1 变化)。
- **单位/坐标**：`Preset.m = 2H/ωs`；发电机为**经典模型**（`E'∠δ` 恒定，`flag_xd=0` 时 epu=Vg）；
  角度多在 **COI 坐标**；`d = m * damping_ratio`。
- **注释/输出**：v1 与 v2 的注释=中文操作说明 + **半角标点**；v1 模块的程序输出=英文；cli/config
  的中文 UI 可用中文。可运行脚本开头用 `SetConsoleOutputCP(65001)` 修 Windows 中文乱码。
- **3 机专属**：v1 的 `reduced_gradient`、`find_reduced_equilibria`、稳定域(2D/3D)、相平面、能量分解
  **假设 ngen=3**；`matlab_xval.reconstruct_cuep` 也是 3 机专用。39 母线现支持**静态初始化 + 时域精确 CCT
  + 能量法 LEA CCT**（后者靠 `bcu_v2/cuep.py` 的通用 closest-UEP 法，见第 8 节方向1）。稳定域/相平面/能量
  分解可视化仍仅 3 机。
- **验证优先级**：物理不变量/解析金标准（P0）> 独立工具（scipy/MATPOWER）> **与 MATLAB 一致（T3）**。
  “与 MATLAB 一致”只代表与该实现一致，不单独等于“正确”。
- **MATLAB 平台只读**：用户明确要求，别写入 `matlab_platform/`。T3 只加载其导出的 `.mat`。

## 7. 已知问题 / 坑（v1 缺陷 + v2 现状）

### v1 六个隐患（2026-09-01 已全部就地修复到 v1；`test_fixes.py` 改为正确性回归 5/5）
> 修复前 v2 靠 `fixes.py` / `cuep.py` 绕开; 现已回灌到 v1 本身。修复后**副效果**: v2 cli 的 mode 段
> (原"噪声陷阱")也自动变正确(9 母线 LEA 0.2274 / 39 母线 0.0847, 退出点均 199)。
1. `solve_sep` 副作用 → 加 `inplace` 开关(默认 True 保留 `build_static_result` 回填); `solve_cuep_from_guess`
   用 `inplace=False`。
2. `run_bcu_experiment` 污染 `postfault.SEP` → 因 1 已解决(CUEP 求解 inplace=False), 污染=0。
3. `find_exitpoint` 初始伪过零(index=0) → 排除初始伪过零(index 0→199)。
4. `trajectory_stable` "回到 SEP" 判据轻阻尼失效 → 加 `criterion` 参数, 默认 `"bounded"`(有界性)。
5. `spm.solve_algebraic` 冷启动偶发不收敛 → 改 `scipy.root`(hybr)+多初值+牛顿回退(冷启动 5/6);
   余个别病态态需 v2 `spm_dae` 连续法根治(7/7)。
6. `find_mgp` 求错(LEA 偏低) → 新增 `bcu_3m9b/cuep.py`(closest-UEP), `run_bcu_experiment` 用
   `controlling_uep` 替 `find_mgp`, LEA 0.057→0.2274(对 MATLAB 0.2275)。

### 其他边界
- **39 母线 H/Xd' 是示例值**（`systems.NE39_H/NE39_XD1`，接近标准 New England）。定量 CCT 前
  **务必用权威参考核对替换**（如 Pai 1989 附录）。装配机理与 SEP/潮流不依赖其精确性。
- **SPM 病态故障网络**（删连接母线）本身**欠定**，属建模边界，非求解器问题（scipy 也救不了）。
- ~~能量法 CUEP 重构暂仅 3 机~~ → **已通用化到任意 ngen**（`bcu_v2/cuep.py`，closest-UEP 法，39 母线可出 LEA）。
  仍待深化：严格 BCU controlling UEP（exit-point→gradient system）与 shadowing。
- one-axis / GFM 是**单机级正确参考模型**，**尚未**接入多机 BCU/能量函数流水线。

## 8. 如何继续开发（下一步建议，按价值排序）

1. ~~**能量法 CUEP 通用化到多机**~~ ✅ **已完成（2026-08-31, `bcu_v2/cuep.py`）**：closest-UEP 法
   （结构化 MOD 初值 + scipy.root 求梯度系统零点 → 取离 SEP 最近 type-1 UEP）。3 机与 MATLAB CUEP 对到
   1e-11、与 3 机网格法对到 3.8e-14；**39 母线首次出能量法 LEA=0.0847s ≤ REA=0.1224s**。已接入 cli 自检。
   `test_cuep.py` 6/6。**关键发现**：MATLAB 的“CUEP”实为 closest-UEP，非 BCU controlling UEP（exit-point
   法在默认故障下判据退化，见 `cuep.py` 模块注释）。**后续**：严格 BCU（exit-point→gradient system→
   controlling UEP）+ shadowing 仍可作为方法学深化（P3.1）。
2. **参数敏感性 / GFL 对比研究脚本**：用 `cli.py` + `--set` 批量扫参（惯量/阻尼/故障位置），出
   CCT 对比表——这是最容易出论文结果的方向。
3. **one-axis 多机化 + 时域仿真**：把 `models.py` 的单机模型接入多机网络（需网络接口 Id/Iq），
   再谈能量函数扩展。
4. **真 IBR 控制**：GFL（PLL+电流环）/ GFM（VSG）完整模型——独立大工程。
5. **文献基准表**：把 9 母线已发表 CCT 硬编码进 `invariants.py` 做 T2 对照（需要你提供数字）。

**开发流程建议**：每次改完，先 `python run_validation.py && python run_matlab_xval.py &&
python test_config.py`，确保“改进不破坏已验证行为”，再提交。

## 9. 验证结果汇总（当前基线，作回归对照）

| 脚本 | 结果 | 关键数字 |
|---|---|---|
| `run_validation.py` | 8/8 | **SMIB 数值 vs 等面积闭式 CCT 差 1e-6**；潮流残差 9e-15；SEP 2.6e-11；LEA≤REA |
| `run_matlab_xval.py`（T3） | 8/8 | Yred 1e-10；SEP 1e-11；**CUEP 1e-11**；LEA 0.2274 vs 0.2275；REA 0.2433 vs 0.2434 |
| `test_fixes.py` | 5/5 | 五隐患修复逐项 |
| `test_spm_dae.py` | 3/3 | DAE 内部一致 9e-9；**稳健性 DAE 7/7 vs v1 3/7** |
| `test_solvers.py` | 3/3 | scipy SEP vs v1 1e-10，更快 |
| `test_p2.py` | 6/6 | 复现 9 母线；**39 母线 SEP 1.7e-15**；ZIP 极限；one-axis→经典 3.8e-10；GFM P=Pset |
| `test_config.py` | 5/5 | 配置加载/校验/9↔39 母线/含 T3 回归 |
| `test_cuep.py` | 6/6 | **通用能量法 CUEP**：9 母线 vs MATLAB 1e-11、vs 3机网格法 3.8e-14；**39 母线 LEA 0.0847≤REA 0.1224** |

## 10. 仓库信息

| 库 | 内容 | 可见性 | 分支 |
|---|---|---|---|
| `github.com/jianzhantianyake-boop/python_platform_lab` | v1（python_bcu） | PRIVATE | main |
| `github.com/jianzhantianyake-boop/python_platform_lab_v2` | v2（python_bcu_v2） | PRIVATE | main |

- 提交作者为本人，**无第三方/AI 署名**。
- 常规更新：在对应目录 `git add -A && git commit -m "..." && git push`。
- 建库/查可见性用 `gh`（已登录）：`gh repo view <owner>/<repo> --json visibility`。

## 11. 术语速查

- **SEP**：稳定平衡点；**UEP**：不稳定平衡点；**CUEP**：受控 UEP（能量法用它算临界能量）。
- **MGP**：最小梯度点（近似 CUEP 的一种追踪法）。
- **CCT**：临界切除时间；**LEA**：能量法（Lyapunov/能量）估计的 CCT；**REA**：时域真值 CCT。
- **COI**：惯性中心坐标；**SPM**：结构保持模型（负荷母线保留为代数节点的 DAE）。
- **GFL/GFM**：跟网型/构网型变流器。**等面积准则**：SMIB 求临界切除角的闭式方法（本项目金标准）。

---

## 12. 全量交叉验证（2026-09-01，MATLAB ↔ Python v2 逐路径）

> 定位：v2 与 MATLAB 是**两个独立平台**，仿真结果应一致（允许平台差异的微小误差）。交叉验证**不以
> MATLAB 为绝对真理**，只作参考；发现分歧时两边都可能有问题。此前 T3（`run_matlab_xval.py`）只覆盖
> `reduced_cct` 一条路径，本轮把覆盖扩到 5/8 路径。详见 `验证覆盖矩阵_CN.md`。

**基础设施**：`matlab -batch`（R2024a 在 PATH）+ `matlab_platform/verify/export_*.m`（不改 B3_MM，
`run` 原 mode 抑制画图 + 精简 save 到 `baseline_*.mat`，`.mat` 已 gitignore、可重生）。
Python 侧统一入口 `python run_full_xval.py`。

**已交叉验证 5/8 路径**（`run_full_xval.py` 5/5 + `reduced_cct` 见 T3 8/8）：

| 路径 | 结果 | 备注 |
|---|---|---|
| reduced_cct | T3 8/8 | 1e-10~1e-11 |
| reduced_region | 4.5e-11 | 平衡点集（SEP+type-1 UEP） |
| reduced_numerical | 1.2e-11 | 故障段末端 thetac |
| two_machine_3d | 6.2e-13 | **交叉验证揪出并修复 D2 双重赋值 bug**（源码 line31=0.45/line215=0.5，EP 用 0.5；用户确认以 0.5 为准，改 `experiments._params_two_machine` d2=0.5） |
| two_machine_gfl | 1.4e-11 | 平衡点集 |
| **spm_cct*** | CCT 0.2053 精确复现 | 见下，E_crit 暂用 MATLAB 参考 |

**SPM 能量法（`bcu_v2/spm_energy.py`，本轮新建）**：物理与公式**全部验证**——关键发现 MATLAB SPM 用
**恒阻抗负荷**（并入 `Yfull_mod`、Sload P/Q=0、网络母线 P=Q=0）。已验证：网络代数解 9e-15、发电机功率
SEP 平衡 1.9e-11、5 项 SPM 势能 E_crit=3.3757 误差 0、CCT 机制给定 E_crit 精确复现 0.2053。
**唯一未闭环 = 自足求 SPM CUEP 网络态的分支选择**（SPM 网络方程在 CUEP 处 11+ 解；与 v1 find_mgp 同源
的数值难题）。这是**明确的下一里程碑**，接手者见 `验证覆盖矩阵_CN.md` 的 SPM 小节（含物理约束
E_crit<故障能量峰值、MATLAB 的 MGP 机制路径）。

**待做（8 路径全覆盖）**：SPM 自足 CUEP（#4 收尾）→ spm_numerical/spm_region（#5/#6，依赖 #4）。

## 13. 交接给 Codex（接手清单）

- **先跑**：`python run_validation.py && python run_matlab_xval.py && python test_config.py && python
  test_cuep.py && python run_full_xval.py`（前四个是回归基线，最后是全量交叉验证 5/5）。
- **本轮新增/改动文件**：
  - v2：`bcu_v2/spm_energy.py`（新，SPM 能量法）、`run_full_xval.py`（新，全量交叉验证）、
    `验证覆盖矩阵_CN.md`（新，覆盖现状+SPM 下一步）。
  - v1：`bcu_3m9b/experiments.py`（`_params_two_machine` d2=0.45→0.5，交叉验证修 D2 bug）。
  - matlab_platform：`verify/export_{region,twomachine,numerical,spm}.m`（新，参考导出脚本）。
- **下一里程碑（SPM 自足 CUEP）**：目标让 `spm_energy` 独立算出 E_critical（≈3.3757）、不依赖 MATLAB。
  两条路：① 移植 MATLAB 的 `Fun_Cal_MGP_SPM`+`Fun_AEiteration_SPM`（DAE 梯度系统沿物理轨迹播种
  `Fun_SEPfslove_SPM` 联合平衡）；② 设计 SPM 版 closest-UEP + 分支选择准则（约束 E_crit<故障能量峰值，
  且 CUEP 网络态与 SEP 连续）。公式/网络解/CCT 机制已全部就绪，只差此环。
- **重跑生成 MATLAB 参考**：`matlab -batch "run('matlab_platform/verify/export_<path>.m')"`。

---

*本文档随代码演进更新。改动核心方程/约定时，请同步更新第 6、7、9 节。*
