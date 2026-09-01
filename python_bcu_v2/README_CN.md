# BCU 平台 v2 · 验证与数值增强

在冻结的 v1(`../python_bcu/bcu_3m9b`)之上,完成 **P0(正确性/可信度)** 与
**P1.2(事件驱动精确 CCT)**。不改 v1 一行,通过 `import bcu_3m9b` 复用。

## 目录结构

```
python_bcu_v2/
  bcu_v2/
    __init__.py      # 路径引导: 把 ../python_bcu 加入 sys.path
    smib.py          # SMIB + 等面积准则闭式 CCT(验证金标准)
    cct.py           # P1.2: 事件驱动 + 二分的精确 CCT(SMIB / 3 机约简)
    fixes.py         # P0.2: 5 个隐患的修正版包装(不改 v1)
    invariants.py    # P0.1: T1 物理不变量 + T2 独立参照 检查集
    spm_dae.py       # P1.3: 严格 DAE 级 SPM(约束流形降阶 + 连续法)
    solvers.py       # P1.1: SciPy 求解器统一层(向后兼容开关)
    systems.py       # P2.1: 通用装配层 + 39 母线 10 机动态案例
    loads.py         # P2.3: ZIP 负荷模型
    models.py        # P2.2 one-axis 同步机 / P2.4 GFM 下垂(单机级参考模型)
    matlab_xval.py   # T3: 与已验证 MATLAB 平台逐层交叉验证(只读参考)
    config.py        # 集中配置: 加载 config.yaml / 校验 / 建 static / 摘要 / 快照
  config.yaml        # 【唯一需要编辑的参数文件】(镜像 MATLAB bcu_config.m)
  cli.py             # 命令行入口(镜像 MATLAB run_bcu.m): run/show/list/validate/xval
  run_validation.py  # 一键跑验证套件 + P1.2 精确 CCT 对照
  run_matlab_xval.py # T3 一键交叉验证
  test_fixes.py test_spm_dae.py test_solvers.py test_p2.py test_config.py  # 各阶段验证
  改进与提升路线图_CN.md
```

## 命令行操作入口(推荐)

```bash
python cli.py show          # 打印当前配置摘要
python cli.py list          # 列出可选 mode / case
python cli.py run           # 用 config.yaml 跑实验(含运行后自检)
python cli.py validate      # P0 物理不变量 + 等面积金标准
python cli.py xval          # T3 与 MATLAB 逐层交叉验证
```
**调参**: 编辑 `config.yaml`(有详细中文注释), 或命令行临时覆盖:
```bash
python cli.py run --mode reduced_region --grid 15
python cli.py run --case case39_modified --auto-params   # 一键切 39 母线(faultline 自动切 [16,17])
python cli.py run --fault-line 8,9 --tunit 5e-4
python cli.py run --set Pm=[0.9,1.3,0.95]                 # 通用 key=value 覆盖
```

## 各阶段独立验证脚本

```bash
python run_validation.py   # P0.1 验证套件(8 项)+ P1.2 精确 CCT
python run_matlab_xval.py  # T3 与 MATLAB 逐层交叉验证
python test_fixes.py       # P0.2 五个修复
python test_spm_dae.py     # P1.3 严格 DAE 级 SPM
python test_solvers.py     # P1.1 scipy 求解器
python test_p2.py          # P2 建模扩展
python test_config.py      # 配置系统 + T3 回归
```

需要 numpy + scipy(已装)。

## 已完成 · 结果

### P0.1 验证套件(zero-MATLAB,8/8 通过)
| 检查 | 结果 |
|---|---|
| 潮流功率失配 → 0 | 9.3e-15 |
| SEP 是真平衡点(残差→0) | 2.6e-11 |
| SEP 是稳定结点(Re λ<0) | max Re λ = -1.34 |
| 总能量单调不增(阻尼耗散) | 净变化 -5.3e-2 |
| 极限:零故障稳 & 长故障失稳 | ✓ |
| 轨迹:v1 RK4 vs scipy solve_ivp | 差 8.1e-4 |
| CCT 夹逼:LEA ≤ 精确 REA | 0.199 ≤ 0.243 ✓ |
| **SMIB:数值 CCT vs 等面积闭式** | **误差 1.0e-6**(金标准) |

> 最硬的一条:SMIB 数值精确 CCT `0.13734s` 与等面积**闭式解** `0.13734s` 吻合到 1e-6,
> 用不依赖任何实现的解析解证明了 CCT 算法本身正确。

### P1.2 事件驱动的精确 CCT
| 方法 | 3 机约简 CCT |
|---|---|
| 能量法 LEA | 0.199 s |
| v1 网格 REA(21 点) | 0.000 s(已知 bug) |
| **P1.2 事件驱动精确** | **0.243 s**(二分 + solve_ivp 事件) |

比固定网格精确 1~2 个数量级,可作为论文级 CCT 数值。

### P0.2 五个隐患修复(5/5 通过,对比 v1 原行为)
1. `solve_sep_pure` — 无副作用(状态精确不变)。
2. `run_experiment_clean` — 零污染(v1 会改写 postfault.SEP)。
3. `find_exitpoint_fixed` — 排除初始伪过零(v1 index=0 → fixed index=199)。
4. `is_stable_bounded` — 轻阻尼下的有界性判据。
5. `solve_algebraic_robust` — scipy.root。**诚实说明**:热启动下与 v1 牛顿成功率相同
   (5/6),单纯换求解器不是万能药;真正的稳健化由 P1.3 的连续法完成(见下)。

### P1.3 严格 DAE 级 SPM(连续法,3/3 通过)
把结构保持模型当"约束流形上的降阶 index-1 DAE":scipy 自适应积分器推进发电机,每次 RHS
求值用**连续法**(上一步解热启动 + scipy.root 校正 + 回退)把负荷代数方程解到机器精度。
| 检查 | 结果 |
|---|---|
| DAE 内部一致(RK45 vs Radau 两积分器) | 差 **9e-9**(积分器无关,解正确) |
| **稳健性(小故障应全稳)** | **DAE 7/7 vs v1 3/7** —— 根治 v1 数值失败 |
| 刚性 Radau 积分可用 | ✓ |

> 关键:v1 在 tc=0.025/0.03/0.08/0.12s 等清除态因代数解冷启动失败, 连续法**全部**解决。

### P1.1 SciPy 求解器统一层(3/3 通过)
`solvers.py`:`nlsolve`(scipy.root / v1 牛顿可开关)、`solve_sep_scipy`(无副作用)、`integrate`
(solve_ivp)。scipy SEP 与 v1 一致到 **1e-10**,且更快(1.0ms vs 2.7ms),坏初值下同样稳健。

### P2.1 通用装配层 + 39 母线(2/2 通过)
`systems.py`:`build_preset` 把"每机 H、Xd'"按 case 发电机顺序装成 v1 Preset(m=2H/ωs),让
`build_static_result` 能跑**任意** case。
- `build_preset` 精确复现 v1 的 9 母线默认参数;
- **`case39_dynamic()` → 39 母线 10 机静态初始化跑通**(Yred 10×10,post SEP 残差 1.7e-15)。
- 平台由此从 9 母线扩展到 **39 母线**。⚠️ 39 母线 H/Xd' 为**示例值**,定量 CCT 前请核对替换。

### P2.3 ZIP 负荷(2/2 通过)
`loads.py`:`zip_load_power`、`zip_algebraic_residual`。三种极限(Z/I/P)正确;aP=1 时残差
**精确等于** v1 的恒功率残差。

### P2.2 / P2.4 单机级参考模型(2/2 通过,基础,未接入 BCU 流水线)
`models.py`:
- **one-axis(磁链衰减)同步机**:在 Xd=Xd' 极限下**精确退化为经典 SMIB**(δ 差 3.8e-10);
- **GFM 下垂逆变器**:稳态 P=Pset 精确。
> 诚实说明:这两个是**单机无穷大母线级**的正确参考实现,已在极限/稳态下校验;**尚未**接入
> 多机 BCU/CCT/稳定域流水线与能量函数——那是更大的后续工作(见路线图 P2.2/P2.4)。

### T3 与 MATLAB 平台逐层交叉验证(8/8 通过)
`matlab_xval.py` 读取已验证 MATLAB 平台导出的参考 `matlab_platform/verify/baseline_reduced.mat`
(只读, 不运行/不改 MATLAB), 与 Python 逐层比对:
| 层 | Python vs MATLAB |
|---|---|
| Yred 预/故障后 | 1.5e-10 / 1.9e-10 |
| SEP 预/故障后角度 + 速度 | 3.3e-11 / 4.4e-11 / 8.7e-11 |
| **CUEP(v2 重构)** | **1.0e-11** |
| **LEA CCT(v2)** | py 0.2274s vs matlab 0.2275s |
| **REA CCT(P1.2)** | py 0.2433s vs matlab 0.2434s |

> 关键发现:核心层(潮流/Yred/SEP)与 MATLAB 吻合到 **1e-10~1e-11**;REA(P1.2)与 MATLAB 吻合
> 到 ~1e-4。**v1 原本的 LEA 偏低是因为 find_mgp 求错**(落在 SEP 附近);v2 用 `reduced_region`
> 的 type-1 UEP 重构 CUEP 后,LEA 也与 MATLAB 对上——T3 既证明了核心正确,又定位并修正了该缺陷。
> 诚实说明:与 MATLAB 一致只代表"与该实现一致";结合 P0 物理不变量 + 等面积金标准才是完整可信度。

### 集中配置 + 命令行入口
`config.yaml`(唯一编辑入口, 详细中文注释)+ `cli.py`(run/show/list/validate/xval)+ `config.py`
(加载/校验/建 static/摘要/快照)。镜像 MATLAB 的 bcu_config.m / run_bcu.m UX:
- 一处改参数(m/damping/Pm/xd1/E、故障、Tfault/Tunit、网格…),或命令行临时覆盖;
- **一键切 9↔39 母线**(`--case case39_modified --auto-params`,39 母线 10 机 SEP 残差 1.7e-15);
- 运行后自检(SEP 残差 + 能量法保守性 LEA≤REA)。
- 无 PyYAML 时内置极简 YAML 解析, 免额外依赖。

## 未做(需要额外条件)
- 文献基准表(9 母线已发表 CCT)需你提供数字。
- 3 机专属可视化(2D/3D 稳定域、能量分解、相平面)在 39 母线不适用;39 母线现支持静态初始化 +
  时域精确 CCT + **能量法 LEA CCT**(2026-08-31, `bcu_v2/cuep.py` 通用 closest-UEP 法, 见下)。
- SPM 病态**故障网络**(删连接母线)本身欠定, 属建模边界而非求解器问题。

## 与 v1 的关系
v1 只读冻结、作为回归基线;v2 只增量、每步用 `run_validation.py` 把关。要把某个 v2 修复
"回灌"进 v1,需你确认后再解冻 v1。
