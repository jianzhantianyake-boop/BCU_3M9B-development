# BCU 3M9B Python 平台 · 操作指南

本指南介绍如何使用当前的 Python 仿真平台（`python_bcu/`），以及在它上面可以
开展哪些研究与实验。核心不依赖 MATLAB；数值求解已引入 SciPy（稳健求解器：SPM 负荷
代数解、CUEP 求根），绘图用 Matplotlib（按需导入）。

> 2026-09-01 更新：v1 已解除冻结，6 个已知隐患全部就地修复，并新增**配置驱动入口**
> `cli.py` + `config.yaml`（见 0.6）。能量法 CUEP 已由 closest-UEP（`bcu_3m9b/cuep.py`）
> 修正，LEA CCT 与 MATLAB 平台一致到 ~1e-11。

> 定位：这是把原 MATLAB 工程（`B3_MM/`）改写成显式 Python 数据结构的**教学/研究
> 实验平台**。核心链路（潮流 → 网络约简 → SEP → 暂态轨迹 → 能量法/时域 CCT →
> 结构保持模型）已可端到端运行；**核心链路已通过 v2 的 T3 与 MATLAB 平台逐层交叉验证到
> 1e-10~1e-11**（v2 复用 v1 这套计算代码，见文末第 7 章）。

---

## 0. 最快上手：中央控制台 main.py（推荐）

不想记 import、也不想分清脚本还是模块？用中央控制台一个入口跑通所有实验：

```bash
cd python_bcu
python main.py            # 弹出中文菜单，输入编号即可
python main.py 3          # 直接跑模式 3（reduced_region）后退出
python main.py a          # 依次跑模式 1-6
python main.py help       # 打印详细操作说明
```

菜单分两部分：8 个**实验模式**（对应 MATLAB 的 `EXPERIMENT_MODE`）+ 若干**工具项**
（`p` 交流潮流 / `s` 静态初始化 / `t` 冒烟自检）。所有模式的图都存到
`python_bcu/figures/`（文件名即模式名）。

### 0.1 八个实验模式（对应 MATLAB EXPERIMENT_MODE）

| 编号 | 模式名 | 对应 MATLAB | 内容 | 输出图 | 状态 |
|---|---|---|---|---|---|
| 1 | `reduced_cct` | Cal_MM_CCT.m | 网络约简：初始化 → 退出点/MGP/CUEP/CCT | δ2-δ3 相平面 | ✅ |
| 2 | `reduced_numerical` | NumSim_MM_Gridframe.m | 网络约简：三段数值仿真（pre/fault/post） | 2×2：角度/速度/相平面/能量 | ✅ |
| 3 | `reduced_region` | Statable_Region.m | 网络约简：二维稳定域（平衡点+分界线） | 稳定域图 | ✅ |
| 4 | `spm_cct` | Cal_MM_CCT_SPM.m | 结构保持：时域 CCT | —（打印） | ✅ |
| 5 | `spm_numerical` | NumSim_MM_Gridframe_SPM.m | 结构保持：数值仿真 | 3 图：角/电压/相平面 | ✅ |
| 6 | `spm_region` | Statable_Region_SPM.m | 结构保持：二维稳定域 | 网格分类+解析分界线 | ✅ |
| 7 | `two_machine_region_3d` | Statable_Region_3D.m | 两机完整模型：三维稳定域 + 稳定流形曲面 | 3D 稳定域图 | ✅ |
| 8 | `two_machine_region_3d_gfl` | Statable_Region_3D_GFL.m | 同上，GFL 参数（阻抗含电阻，低惯量） | 3D 稳定域图 | ✅ |

> 也可以在 Python 交互环境里直接调用：`from bcu_3m9b import experiments as e; e.mode_reduced_region()`。

### 0.2 每个模式细看

**1 `reduced_cct`** — 从预故障 SEP 沿故障网络积分找退出点，用 **closest-UEP**（`cuep.py`，修正
原 `find_mgp` 求错）求 CUEP，给出能量法 LEA CCT 与时域 REA CCT；画 δ2-δ3 相平面（故障轨迹 +
SEP + 退出点 + CUEP）。参考量级：LEA ≈ 0.2274 s，REA ≈ 0.24 s。

```python
from bcu_3m9b import experiments as e
e.mode_reduced_cct()      # 图：figures/reduced_cct_phaseplane.png
```

**2 `reduced_numerical`** — prefault/fault/postfault 三段定步长积分，输出 2×2 多图：
COI 角度、COI 相对速度、δ2-δ3 相平面、故障后能量（势能/动能/总能量）。

```python
e.mode_reduced_numerical(fault_time=0.2, postfault_time=3.0)
```

**3 `reduced_region`** — 忠实移植 MATLAB `f_reducedstate` 梯度系统：在 δ2-δ3 平面网格
搜索平衡点，按雅可比非负特征值个数分类（蓝=稳定 SEP，红=type-1 UEP），对每个
type-1 UEP 反向积分稳定流形画出**稳定域分界线**（separatrix），并周期平铺。这是本组
里最“漂亮”的图，对应经典 BCU 稳定域示意。

```python
e.mode_reduced_region(grid_points=21)   # 网格越大越精细也越慢
```

**4 `spm_cct`** — 结构保持模型的时域 CCT：在 `[0, fault_max]` 扫切除时刻，故障段用约简
模型积分（数值稳健），故障后段用 SPM 正向积分并按“有界性（是否失步）”判稳。

> 注：故障网络删了连接母线，其 SPM 代数方程病态；SPM 代数解偶发数值失败处会回退到
> 约简模型的有界性判据。故这是**可运行近似**，非原版反向 DAE（`ode15s`）。

**5 `spm_numerical`** — 故障段约简 + 故障后 SPM，输出 3 图：发电机 COI 角、负荷母线
电压、δ2-δ3 相平面。

**6 `spm_region`** — 在 δ2-δ3 平面撒初值，正向 SPM 积分后按有界性分类稳定/不稳定，
并叠加约简模型的解析分界线作参考（SPM 稳定点应落在分界线内，两法自洽）。

```python
e.mode_spm_region(grid_points=15)       # 较慢（每个网格点一次 SPM 仿真）
```

**7 `two_machine_region_3d` / 8 `two_machine_region_3d_gfl`** — 独立于 9 母线电网，用
**两机完整模型** `f_2m`（状态 `[δ12, ω12, ω_sum]`）。在 δ12 网格找平衡点、按非负特征值
个数分类，对每个 type-1 UEP 扫其**二维稳定子空间**反向积分出**稳定流形曲面**（即 3D
稳定域边界），type-2 UEP 出一维曲线，并叠加一条“故障 → 切除”样例轨迹。GFL 版只是换一
组参数（阻抗含电阻、低惯量），**不是另一套动力学**（这与 MATLAB `Statable_Region_3D` 和
`_GFL` 两脚本结构相同、仅参数不同一致）。

```python
e.mode_two_machine_region_3d()       # 图：figures/two_machine_region_3d.png
e.mode_two_machine_region_3d_gfl()   # 图：figures/two_machine_region_gfl.png
```

### 0.3 改参数

- 稳定域网格密度：改 `main.py` 顶部 `CONFIG` 的 `region_grid`、`spm_region_grid`。
- 更细的时长/步长：改 `bcu_3m9b/experiments.py` 各 `mode_*` 函数签名的默认值。

### 0.4 与本指南其余章节的关系

下面第 1–8 章讲的是**底层 API**（可在 Python 交互环境里逐个函数调用），本节的
`main.py` / `experiments.py` 则是把这些 API 编排成“一键模式”的上层入口。想快速出结果
看本节；想自定义实验看下面的 API。

### 0.5 已知数值提示

- 本系统阻尼很轻（`d ≈ 0.1·m`），暂态要振荡几十秒才 settle，故稳定性判据用**有界性**
  （轨迹是否失步）而非“末端回到 SEP”。`trajectory_stable` 默认 `criterion="bounded"`。
  参考量级：reduced REA CCT ≈ 0.24 s，spm_cct ≈ 0.225 s；能量法 LEA CCT ≈ 0.2274 s。
- （2026-09-01 起）`run_bcu_experiment` 求 CUEP 用 `inplace=False`，**不再改写** `postfault`
  的 SEP 字段；跨模式共用同一个 `static` 已安全（原“各自 build 以防污染”的约束解除）。

### 0.6 配置驱动入口 cli.py（推荐用于复现 / 批量扫参）

除交互菜单 `main.py`，平台还提供**配置驱动**入口 `cli.py`（镜像 MATLAB `run_bcu.m` / v2 的
`cli.py`），把“改一处参数 → 一键跑 → 自动自检”串起来，适合复现与批量扫参。

```bash
cd python_bcu
python cli.py show          # 打印当前配置摘要(不跑)
python cli.py list          # 列出可选 mode / case
python cli.py run           # 用 config.yaml 跑 mode(含运行后自检)
python cli.py validate      # 基础正确性自检(潮流/SEP 残差 + 能量法 LEA<=REA)
```

**改参数两种方式**：

- 编辑 `config.yaml`（唯一参数文件，有中文注释）：`mode` / 发电机参数（`m`/`damping_ratio`/`Pm`/
  `xd1`/`E`）/ `faultline` / `Tfault` / `Tunit` / `region_grid` / `spm_region_grid` 等。
- 命令行临时覆盖（不改文件）：

```bash
python cli.py run --mode reduced_region --grid 15
python cli.py run --fault-line 8,9 --tunit 5e-4
python cli.py run --set Pm=[0.9,1.3,0.95]      # 通用 key=value 覆盖(值走 JSON 解析)
```

**运行后自检**（`reduced_cct` / `reduced_numerical`）打印：postfault SEP 残差、能量法 CUEP
（type-1 + `V(CUEP)`）、`LEA / REA CCT` 与保守性 `LEA<=REA`。每次 `run` 会把配置快照存到
`results/`（可追溯；`save_snapshot: false` 关闭）。

> `main.py` 与 `cli.py` 调的是同一套底层功能：**想交互探索用 `main.py`，想复现/批量用 `cli.py`**。
> v1 只支持 `case9_v2`（3 机 9 母线）；39 母线 / 更多模型 / T3 交叉验证在 v2（`../python_bcu_v2/`）。

---

## 1. 环境与运行

| 项 | 说明 |
|---|---|
| Python | 3.12（已验证） |
| 依赖 | `numpy`、`scipy`（均必需；scipy 用于稳健求解器）、`matplotlib`（仅绘图时需要） |
| 位置 | `.../ntu_cladue_only/代码/BCU_3M9B-main/python_bcu/` |

**最快的运行方式**（在 `python_bcu/` 目录下）：

```bash
cd python_bcu
python tests/smoke_test.py         # 冒烟自检，应打印 smoke_test: PASS
python scripts/run_static.py       # 9 母线静态初始化，打印电压/SEP/残差
python scripts/run_experiment.py   # 完整 BCU/能量/CCT/暂态实验
python scripts/run_spm.py          # 结构保持模型最小实验
python scripts/run_two_machine.py  # 两机模型 + 平衡点扫描
```

脚本内部会把包根目录加入 `sys.path`，因此**从 `python_bcu/` 目录**直接跑即可。
若想在任意目录 `import bcu_3m9b`，在 `python_bcu/` 下执行 `pip install -e .`（有
`pyproject.toml`）。

---

## 2. 模块地图（`bcu_3m9b/`）

| 模块 | 职责 | 关键函数/类 |
|---|---|---|
| `types.py` | 全部数据结构（替代 MATLAB base workspace） | `CaseData` `PFData` `Preset` `BaseValue` `NetworkState` `Trajectory` `StaticResult` `CCTResult` |
| `cases.py` | 内置案例数据 | `case9_v2()` `case39_modified()` |
| `powerflow.py` | 交流潮流 + pfdata 转换 | `solve_power_flow` `to_pfdata` `generator_internal_emf` |
| `network.py` | 导纳组装、故障设置、Kron 约简 | `rxb_to_yfull` `add_load_admittance` `kron_reduce` `remove_fault_bus` `remove_fault_line` |
| `equilibrium.py` | SEP/CUEP 残差与求解 | `electrical_power` `solve_sep` `sep_check` `solve_cuep_from_guess` |
| `dynamics.py` | 暂态轨迹、退出点、时域 CCT | `integrate_reduced` `find_exitpoint` `trajectory_stable` `time_domain_cct` |
| `energy.py` | 势能/动能/MGP/能量法 CCT | `potential_energy` `trajectory_energy` `energy_cct` `find_mgp` |
| `cuep.py` | 通用能量法 CUEP + LEA（closest-UEP，任意 ngen） | `coi_mismatch` `find_type1_ueps` `controlling_uep` `energy_lea_cct` |
| `spm.py` | 结构保持模型（显式 DAE 近似；scipy 稳健代数解） | `solve_algebraic` `simulate_spm` |
| `two_machine.py` | 两机模型 + 平衡点分类 | `TwoMachineParameters` `f_2m` `equilibria` `simulate_two_machine` |
| `stability_region.py` | 稳定域/向量场批量实验 | `vectorfield_norm` `simulate_grid` |
| `numerics.py` | 数值内核（无 SciPy） | `newton_solve` `numerical_jacobian` `rk4_step` |
| `bcu.py` | 一键入口 | `default_preset` `build_static_result` `run_bcu_experiment` |
| `config.py` | 集中配置（`cli.py` 用） | `load_config` `validate_config` `build_static_from_config` `save_snapshot` |
| `matlab_compat.py` | MATLAB 同名函数薄封装（`Fun_*`） | 便于与原文件逐个对照 |
| `plotting.py` | 可选绘图 | `plot_trajectory` `plot_energy` |

顶层入口（可从包直接导入）：

```python
from bcu_3m9b import (
    case9_v2, case39_modified,
    solve_power_flow, to_pfdata,
    default_preset, build_static_result, run_bcu_experiment,
)
```

---

## 3. 核心工作流

### 3.1 一次静态初始化（对应 MATLAB `Cal_MM_Static`）

```python
from bcu_3m9b import build_static_result

static = build_static_result()          # 默认 case9_v2 + default_preset
print(static.prefault.sep_delta)        # 预故障 SEP 功角 (rad)
print(static.postfault.sep_delta)       # 故障后 SEP 功角
print(static.prefault.yred.shape)       # 约简导纳维度 (3, 3)
```

`build_static_result` 一步完成：潮流 → pfdata → 预故障/故障/故障后三套导纳与
Kron 约简 → 预故障与故障后 SEP，并生成结构保持模型所需的发电机优先重排。返回
`StaticResult`（内含 `preset / basevalue / pfdata / prefault / fault / postfault / emf / case`）。

### 3.2 一键完整实验（对应 BCU 主流程）

```python
from bcu_3m9b import build_static_result, run_bcu_experiment

static = build_static_result()
res = run_bcu_experiment(static,
                         fault_time=0.2,     # 故障持续时长 (s)
                         tunit=1e-3,         # 积分步长 (s)
                         postfault_time=2.0, # 故障后仿真时长 (s)
                         cct_samples=11)     # 时域 CCT 网格点数

res["exit_index"], res["exit_time"]   # 退出点
res["cuep_delta"]                     # CUEP (closest-UEP)
res["cuep_result"]                    # CUEPResult(type-1 特征值 / V(CUEP) / 候选 UEP 数)
res["cuep_source"]                    # CUEP 来源标签(closest-UEP / MGP fallback)
res["mgp"]["theta_mgp"]               # MGP 候选点(诊断/回退用)
res["critical_energy"]                # 临界能量 V(CUEP)
res["lea"].cct, res["lea"].flag_cct   # 能量法 (LEA) CCT
res["rea_cct"]                        # 时域网格 (REA) CCT
```

---

## 4. 可开展的研究与实验

下面每一条都是当前平台**已能跑**的实验方向，并给出最小代码。

### 4.1 交流潮流研究
研究不同案例/负荷/发电设定下的潮流解、母线电压分布、支路潮流。

```python
from bcu_3m9b import case9_v2, case39_modified, solve_power_flow, to_pfdata
pf = solve_power_flow(case39_modified(), tol=1e-8)
print(pf.success, pf.iterations, pf.residual_norm)
pfd = to_pfdata(pf)
print(pfd.voltage)          # [幅值, 角度(deg)]
print(pfd.branch_powerflow) # 支路两端 P/Q
```
可做：修改 `case9_v2()` 的 bus/gen/branch 观察潮流灵敏度；对比 9 母线与 39 母线。

### 4.2 平衡点（SEP）与临界点（CUEP）
研究稳定平衡点位置、随参数的漂移，以及以 MGP 为初值搜索 CUEP。

```python
from bcu_3m9b import build_static_result
from bcu_3m9b.equilibrium import solve_sep, solve_cuep_from_guess
s = build_static_result()
# 故障后 SEP 已在 static 内；也可对任意工况重解
delta, wpu, perr, ok, it = solve_sep(s.preset, s.postfault, s.basevalue)
```
可做：扫描机械功率 `preset.pmpu`、内电势 `preset.epu` 观察 SEP 变化与失稳边界。

### 4.3 暂态轨迹与退出点
研究故障期/故障后的功角摇摆、COI 坐标下的相对运动、退出点判据。

```python
from bcu_3m9b.dynamics import integrate_reduced, find_exitpoint, trajectory_stable
import numpy as np
s = build_static_result()
d0 = s.prefault.sep_delta
w0 = np.full(s.preset.ngen, s.prefault.sep_omegapu * s.basevalue.omega_b)
traj = integrate_reduced(0.2, 1e-3, s.fault, s.preset, s.basevalue, d0, w0)
idx = find_exitpoint(traj, s.postfault, s.preset)
stable = trajectory_stable(traj, s.postfault, s.preset)
```
`integrate_reduced` 的 `semi_rk4=True` 复刻原 MATLAB 步进顺序，`False` 用标准 RK4，
可做**数值方法对比**教学实验。

### 4.4 能量函数与临界能量（LEA）
研究势能三分解、动能、总能量沿轨迹的演化，以及能量法 CCT。

```python
from bcu_3m9b.energy import trajectory_energy, potential_energy, energy_cct
en = trajectory_energy(traj, s.preset, s.postfault)
en["ep"], en["ek"], en["total"]     # 逐步能量
```
可做：比较不同 `preset.path_energy_cal`（0 解析 / -1 忽略电导项 / >0 多段梯形积分）
对势能路径项 Ep3 的影响。

### 4.5 CCT 双路径对比（能量法 LEA vs 时域 REA）
用同一工况同时得到能量法 CCT 与时域网格 CCT，研究两者差异。

```python
res = run_bcu_experiment(s, fault_time=0.2, tunit=1e-3, postfault_time=2.0, cct_samples=21)
print("LEA:", res["lea"].cct, "REA:", res["rea_cct"])
```
可做：加密 `cct_samples`、缩小 `tunit` 观察 REA 收敛；研究保守性差异。

### 4.6 CUEP 求解 / MGP 边界追踪
研究受控不稳定平衡点（CUEP）的定位与临界能量。**推荐用 `cuep` 模块**（closest-UEP，任意 ngen，
修正原 `find_mgp` 求错）：

```python
from bcu_3m9b.cuep import controlling_uep, energy_lea_cct, find_type1_ueps
cres = controlling_uep(s)             # 离 SEP 最近的 type-1 UEP
cres.cuep, cres.v_cuep, cres.eig_reduced   # CUEP 角 / 临界能量 V(CUEP) / 约简特征值
lea = energy_lea_cct(s)               # 能量法 LEA CCT
ueps = find_type1_ueps(s)             # 全部 type-1 UEP(按离 SEP 距离排序)
```

原 `find_mgp` 仍保留，仅作诊断/回退：

```python
from bcu_3m9b.energy import find_mgp
mgp = find_mgp(traj.thetac[idx], s.postfault, s.preset)
mgp["theta_mgp"], mgp["found"], mgp["norm"]
```

### 4.7 结构保持模型（SPM）
研究不消去负荷节点的结构保持网络下发电机 + 负荷母线的联合动态。

```python
from bcu_3m9b.spm import simulate_spm
import numpy as np
w0 = np.full(s.preset.ngen, s.postfault.sep_omegapu * s.basevalue.omega_b)
out = simulate_spm(0.05, 2e-3, s.postfault, s.preset, s.basevalue,
                   s.postfault.sep_delta, w0)
out["delta_coi"][-1], out["algebraic"][-1]   # 末端功角 / 负荷节点代数解
```
> 注意：SPM 用「每步先牛顿求代数解、再 RK4 推进发电机」的可运行近似，不等同于
> MATLAB `ode15s` 的严格 DAE 逐点轨迹。

### 4.8 两机模型与平衡点分类
经典两机等值系统：相平面、平衡点稳定性、故障积分。

```python
from bcu_3m9b.two_machine import TwoMachineParameters, equilibria, simulate_two_machine
```
`equilibria` 返回每个平衡点的雅可比特征值与非负特征值个数（稳定性分类）。

### 4.9 稳定域 / 向量场批量实验
对一批初值批量积分并粗标稳定/发散，勾勒吸引域。

```python
from bcu_3m9b.stability_region import simulate_grid, vectorfield_norm
```

### 4.10 绘图（可选，需 matplotlib）
```python
from bcu_3m9b.plotting import plot_trajectory, plot_energy
plot_trajectory(traj); plot_energy(en)
```

---

## 5. 关键参数与开关（`Preset`）

| 字段 | 含义 | 常用取值 |
|---|---|---|
| `m` / `d` | 惯量 / 阻尼向量 | 见 `default_preset()` |
| `pmpu` / `epu` | 机械功率 / 内电势（pu） | 扫描研究失稳边界 |
| `xd1` / `flag_xd` | 暂态电抗 / 是否反算内电势 | `flag_xd=0` 直接用端电压 |
| `fault_line` / `fault_position` | 故障线路 `[i, j]` / 故障母线索引 | 默认 `[9,6]`, `0`→母线 9 |
| `path_energy_cal` | 势能路径项 Ep3 算法 | `0` 解析 / `-1` 忽略 / `>0` 多段梯形 |

改参数的推荐做法：复制默认再改。

```python
from bcu_3m9b import default_preset, build_static_result
p = default_preset(); p.fault_line[:] = [7, 5]
s = build_static_result(preset=p)
```

---

## 6. 与 MATLAB 对照阅读

`matlab_compat.py` 提供与原 `Fun_*` 同名的薄封装（如 `Fun_RXB2Yfull`、
`Fun_Yfull2Yred`、`Fun_SEPfslove`、`Fun_TrajIter_SRF`、`Fun_Cal_Exitpoint`、
`Fun_Cal_PotentialEnergy`、`Fun_Cal_MGP`、`Fun_Cal_CCT_Energy`），便于把 `B3_MM/`
里的 MATLAB 文件逐个对照到 Python 实现。函数只用显式传参，不读全局变量。

---

## 7. 已知限制与可扩展方向

- **定量交叉验证已建立（在 v2）**：v1 的核心链路（潮流/Yred/SEP/CUEP/CCT）已由 v2 的 T3
  （`bcu_v2/matlab_xval.py`，读 MATLAB 平台导出的 `verify/baseline_reduced.mat`）逐层对齐到
  **1e-10~1e-11**（Yred 1e-10、SEP 1e-11、CUEP 1e-11、LEA 0.2274 vs 0.2275、REA 0.2433 vs 0.2434）。
  因 v2 `import bcu_3m9b`，被验证的就是 v1 这套计算代码；T3 脚本本身在 v2（`../python_bcu_v2/`）。
  v1 潮流用自带有限差分牛顿（非 MATPOWER `runpf`），数值仍与 MATLAB 吻合。
- **SPM 简化**：SPM 三模式为“故障段约简 + 故障后 SPM 正向积分”的可运行近似，非原版
  反向 DAE（`ode15s`）逐点复刻；未做故障切换、ZIP 负荷、严格 DAE 积分器。
- **CUEP 为 closest-UEP**：`run_bcu_experiment` 用 `cuep.controlling_uep`（离 SEP 最近的 type-1 UEP）
  求 CUEP，修正了原 `find_mgp` 求错致 LEA 偏低（LEA 0.057→0.2274，对 MATLAB 0.2275）；仅在
  closest-UEP 搜索失败时才回退 MGP 近似。closest-UEP 是 controlling UEP 的常用工程近似，定量结论
  建议结合时域 REA 交叉核对。
- **绘图已覆盖 8 模式**：第 0 章已提供 δ2-δ3 相平面、三段数值多图、二维稳定域（含解析
  分界线）、两机三维稳定域 + 稳定流形曲面（含 GFL 参数）等；**尚未做**阻尼能量分组、
  三维曲面着色/透明填充等更精细的科研样式。
- **39 母线动态实验**：默认最小入口仍用 9 母线；39 母线动态参数/故障设置待补。

---

## 8. 一分钟上手清单

```bash
cd python_bcu
python tests/smoke_test.py   # 1) 确认环境（应打印 PASS）
python main.py               # 2) 启动中央控制台，输入编号跑模式（见第 0 章）
python main.py a             # 3) 一次性跑完模式 1-8，图存到 figures/
python cli.py run            # 3') 或配置驱动: 改 config.yaml 一键跑+自检（见 0.6）
python cli.py validate       #     基础正确性自检（SEP 残差 + LEA<=REA）
```
```python
# 4) 进入交互式探索（底层 API）
from bcu_3m9b import build_static_result, run_bcu_experiment
s = build_static_result()
res = run_bcu_experiment(s)
# 或直接调用某个模式
from bcu_3m9b import experiments as e
e.mode_reduced_region()
```
