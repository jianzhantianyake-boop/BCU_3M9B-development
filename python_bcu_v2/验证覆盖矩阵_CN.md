# BCU 平台 · 验证覆盖矩阵（对 MATLAB 交叉验证现状）

> 目的：把"MATLAB 平台的哪些路径，Python 已经验证过、哪些没有、各自缺什么"一次性钉清楚，
> 消除"T3 8/8 绿"造成的"全平台已验证"错觉。**这是现状清单，不是完成报告。**
> 最后更新：2026-09-01。单位：角度 rad，时间 s，功率/导纳 pu。

---

## 0. 一句话结论

**目前只有 `reduced_cct` 一条路径真正对 MATLAB 交叉验证过（到 1e-10~1e-11）；其余 7 条路径要么
只验证了静态层、要么完全没对过 MATLAB，其中 SPM 能量法在 Python 里根本没实现。**

"T3 8/8 通过"这个数字的**真实含义**是：*在 reduced 这一条链路上，与 MATLAB 全部一致*——它
**不代表**整个平台都验证过了。

---

## 1. 关键澄清：现有 T3 的真实边界

`bcu_v2/matlab_xval.py`（即 `python cli.py xval`）读的基准文件是
`matlab_platform/verify/baseline_reduced.mat`——**文件名就是 `reduced`**。它的 8 个验证层
**全部属于 reduced（网络约简）模型**：

| T3 的 8 层 | 模型 |
|---|---|
| Yred 预/故障后（2） | reduced |
| SEP 预/故障后角度+速度（3） | reduced |
| CUEP（1） | reduced |
| LEA CCT（1） | reduced 能量法 |
| REA CCT（1） | reduced 时域 |

**没有任何一层触及 SPM 或 two_machine。**

> 另注：`matlab_platform/verify/` 下的 `v0_baseline.m` / `v1_unit.m` / `v2_energy.m` 是 **MATLAB
> 平台自身**的分层自检（基线快照 / 底层函数解析解 / 能量守恒律），属于 **MATLAB 侧内部验证**，
> **不是** Python↔MATLAB 的跨平台交叉验证。别把它们算进"Python 已验证"。

---

## 2. 覆盖矩阵（主表）

MATLAB 8 条路径（`run_bcu.m` 的 `switch cfg.mode`）逐条对照：

> **进度更新 2026-09-01**：阶段 A+B 完成，5/8 路径已对 MATLAB 交叉验证（`run_full_xval.py`）。
> 交叉验证过程揪出并修复了 two_machine 的 D2 双重赋值 bug。SPM 三条经评估为"完整 DAE 管线待移植"。

| # | MATLAB 路径 | MATLAB 脚本 | Python 实现 | MATLAB 参考 | 已交叉验证 | 误差/缺口 |
|---|---|---|---|---|---|---|
| 1 | `reduced_cct` | `Cal_MM_CCT.m` | ✅ 完整 | ✅ `baseline_reduced.mat` | ✅ **T3 8/8** | 1e-10~1e-11 |
| 2 | `reduced_numerical` | +`NumSim_MM_Gridframe.m` | ✅ 完整 | ✅ `baseline_numerical.mat` | ✅ **PASS** | 故障段末端 1.2e-11 |
| 3 | `reduced_region` | +`Statable_Region.m` | ✅ 完整 | ✅ `baseline_region.mat` | ✅ **PASS** | 平衡点集 4.5e-11 |
| 4 | `spm_cct` | `Cal_MM_CCT_SPM.m` | ✅ 能量法(`spm_energy.py`) | ✅ `baseline_spm.mat`(0.2053) | 🟡 **部分**(见下) | 自足 CUEP 分支待做 |
| 5 | `spm_numerical` | +`NumSim_MM_Gridframe_SPM.m` | ⚠️ 近似 | ❌ 无 | ❌ | 依赖 #4 |
| 6 | `spm_region` | +`Statable_Region_SPM.m` | ⚠️ 近似 | ❌ 无 | ❌ | 依赖 #4 |
| 7 | `two_machine_region_3d` | `Statable_Region_3D.m` | ✅ 完整 | ✅ `baseline_twomachine.mat` | ✅ **PASS**(修 D2) | 平衡点集 6.2e-13 |
| 8 | `two_machine_region_3d_gfl` | `Statable_Region_3D_GFL.m` | ✅ 完整 | ✅ `baseline_twomachine_gfl.mat` | ✅ **PASS** | 平衡点集 1.4e-11 |

### SPM（#4）当前状态：物理全验证，自足 CUEP 分支待闭环
`bcu_v2/spm_energy.py` 已逐块对 MATLAB 交叉验证（关键发现：MATLAB SPM 用**恒阻抗负荷**，负荷并入
`Yfull_mod`、`Sload` 的 P/Q=0、网络母线 P=Q=0）：

| SPM 组件 | 验证结果 |
|---|---|
| `solve_spm_network`（网络代数解，恒阻抗） | SEP 网络态 **9e-15** |
| `spm_generator_power`（经 Yfull_mod） | SEP COI 功率失配 **1.9e-11** |
| `spm_potential_energy`（5 项含网络电压） | 用 MATLAB 网络态得 E_crit=3.3757，**误差 0** |
| `spm_fault_energy_cct`（CCT 机制） | 给定 E_crit=3.3757 → **CCT 0.2053 精确复现** |

**唯一未闭环：自足求 SPM CUEP 网络态。** 这是真正的数值难题（非 MATLAB bug，非公式错）：
- SPM 网络方程在 CUEP 发电机角处有 **11+ 个解**（能量 −6 到 +9.5）；要选"控制 UEP"那个物理分支。
- 直线/连续法（从 SEP 或退出点）都跨过分支边界，落到非物理解（E_crit=7.58，> 故障能量峰值 5.568，
  会给出"永远稳定"——物理上不可能，故被否定）。
- reduced 梯度系统（=SPM 梯度系统）从退出点追踪 **不干净收敛到 CUEP**（停在 δ3≈2.89 而非 2.0）——
  这与当年 **v1 `find_mgp` 求错同源**；reduced 路径当初用 **closest-UEP** 绕过。
- MATLAB 靠 `Fun_Cal_MGP_SPM` + `Fun_AEiteration_SPM` 沿物理轨迹连续跟踪播种 `Fun_SEPfslove_SPM`
  联合平衡，锁定分支。

**下一里程碑（交接给接手者）**：移植 MATLAB 的 SPM MGP 机制（DAE 梯度系统 + 网络准代数），
或为 SPM 设计稳健的 closest-UEP 式分支选择准则（约束：E_crit < 故障能量峰值），让 `spm_energy` 自足
求出 E_critical。物理约束与全部公式已就绪，只差这一个分支选择环节。

> 历史注记：曾误判 MATLAB 有 deltacoi bug（`CUEP_net_theta` 残差 11.3）。后用物理约束
> （E_crit 须 < 故障能量峰值 5.568）判定：MATLAB 的 3.3757 物理站得住（保守，LEA 0.2053 < REA 0.2433），
> 我此前算的 7.58 才是选错分支。交叉验证不以 MATLAB 为绝对真理，但此处 MATLAB 结果经独立物理判据成立。

---

## 3. 缺口分三类（决定各自工作量）

- **A 类｜只缺 MATLAB 参考对照**（Python 已完整实现，跑一次 MATLAB 导出关键量即可比）：
  路径 2（轨迹级）、3、7、8。**工作量小**——阶段 A（导出基础设施）+ 阶段 B（对比脚本）。
- **B 类｜缺 Python 实现**（Python 是简化/近似，得先补算法）：
  路径 4（`spm_cct` 的 SPM 能量法四件套：Exitpoint_SPM / MGP_SPM / CUEP_SPM / CCT_Energy_SPM）、
  可能连带 5。**工作量大**——阶段 C。
- **既有的**：路径 1，已闭环，作回归基线。

> SPM 与 reduced 的能量法是**两套独立管线**：SPM 势能含**网络电压项**，CUEP 由 fsolve 同时满足
> 发电机+网络节点平衡。所以 MATLAB `spm_cct` 的 `CCT(LEA)=0.2053s` 与 reduced 的 `0.2275s`
> **本就不同**（不同模型的能量法），不能互相顶替。Python 目前的 `spm_cct` 借用了 reduced 的能量
> 法数字，是错位的。

---

## 4. 各路径当前可信度分级（诚实标注）

| 路径 | 可信度 | 依据 |
|---|---|---|
| `reduced_cct` | **高** | MATLAB 1e-11 + 物理不变量(P0) + 等面积闭式金标准(1e-6) 三重背书 |
| `reduced_numerical` | 中 | 静态初始化经 T3 验证；三段轨迹与 MATLAB 未逐点对比 |
| `reduced_region` | 中 | 平衡点/分界线算法忠实移植 MATLAB `f_reducedstate`，但无数值对照 |
| `spm_cct` | **低** | 简化实现（缺 SPM 能量法）；从未对 MATLAB；时域近似≈0.225s |
| `spm_numerical` / `spm_region` | **低** | 可运行近似（非严格反向 DAE）；从未对 MATLAB |
| `two_machine_region_3d(_gfl)` | 低-中 | 实现完整、逻辑与 MATLAB 同构，但无数值对照 |

> 注：`test_spm_dae.py`（DAE 7/7 vs v1 3/7）是 **zero-MATLAB 的内部自洽性**验证（两积分器一致 +
> 稳健性），**不是**对 MATLAB 的交叉验证——不提升 SPM 路径的"对 MATLAB 可信度"。

---

## 5. 补齐计划（钉现状后按此施工）

前提已确认：**`matlab -batch` 可无界面运行**（R2024a，`setup_bcu_paths` 可用），能生成新参考。

- **阶段 A**：统一 MATLAB `-batch` 导出脚本——每条路径**精简导出**关键数值量（平衡点/CUEP/CCT/
  轨迹末端/网格分类计数等，**不导全轨迹**，避免再出 271MB 巨型 .mat）→ `verify/baseline_<path>.mat`。
- **阶段 B**（A 类，工作量小）：reduced_region、two_machine 3D/GFL、reduced_numerical 交叉验证。
- **阶段 C**（B 类，工作量大）：实现 SPM 能量法四件套 → 对 MATLAB `spm_cct` 交叉验证（目标对上 0.2053s）。
- **阶段 D**：spm_region / spm_numerical 验证。
- **阶段 E**：全量回归脚本 `run_full_xval.py`（8 路径一次跑绿）+ 本矩阵更新为完成状态。

---

## 6. 附：MATLAB 脚本 ↔ Python 实现对照（施工索引）

| MATLAB 关键函数 | Python 对应 | 状态 |
|---|---|---|
| `Cal_MM_Static.m` | `bcu.build_static_result` | ✅ 已验证 |
| `Fun_Cal_Exitpoint.m` | `dynamics.find_exitpoint` | ✅ 已验证 |
| `Fun_Cal_MGP.m` / CUEP | `cuep.controlling_uep`(closest-UEP) | ✅ 已验证(reduced) |
| `Fun_Cal_CCT_Energy.m` | `energy.energy_cct` | ✅ 已验证 |
| `Statable_Region.m` | `experiments.find_reduced_equilibria` + 分界线 | ⚠️ 未对照 |
| `Fun_Cal_Exitpoint_SPM.m` | **无** | ❌ 待实现(阶段C) |
| `Fun_Cal_MGP_SPM.m` | **无** | ❌ 待实现(阶段C) |
| `Fun_Cal_CCT_Energy_SPM.m` | **无** | ❌ 待实现(阶段C) |
| `Statable_Region_3D.m` / `_GFL.m` | `experiments.mode_two_machine_region_3d(_gfl)` | ⚠️ 未对照 |

---

*本矩阵是"全量交叉验证"工程的现状基线。每完成一条路径的交叉验证，回来把对应行的"已交叉验证"
改为 ✅ 并填误差量级。*
