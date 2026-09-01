# BCU 平台 v2 · 当下交互操作手册

> 目的: 回答一个问题——"我现在坐在这个平台前, 应该怎么操作, 每一步会看到什么, 哪些数字能信"。
> 适用范围: `python_bcu_v2/` 当前状态(v1 已解冻并回灌修复, v2 增量)。最后更新: 2026-09-01。
> 单位约定: 角度 rad, 角速度 rad/s, 时间 s, 功率/电压/导纳 pu, 频率 Hz。

---

## 0. 一句话定位

- **唯一日常入口是 `cli.py`**(镜像 MATLAB 的 `run_bcu.m`)。你 90% 的操作就是改 `config.yaml`
  或加命令行参数, 然后 `python cli.py run`。
- **唯一需要编辑的参数文件是 `config.yaml`**(镜像 MATLAB 的 `bcu_config.m`), 有详细中文注释。
- 所有命令都在 `python_bcu_v2/` 目录下执行(v2 靠 `bcu_v2/__init__.py` 把 `../python_bcu`
  加进 `sys.path`, 换目录会 import 失败)。

```bash
cd "C:\Users\WangS\Desktop\26年第二学期学习文件\论文文献\ntu_cladue_only\代码\BCU_3M9B-main\python_bcu_v2"
```

---

## 1. 每次开工先做: 确认平台健康(3 条命令)

改任何东西之前, 先跑这三条确认基线是绿的。三个都过, 说明环境 + 代码 + 参考数据都健康。

```bash
python run_validation.py     # P0: 物理不变量 + 等面积金标准, 期望 8/8
python run_matlab_xval.py    # T3: 与 MATLAB 逐层对比, 期望 8/8
python test_config.py        # 配置系统 + T3 回归, 期望 5/5
```

关键基线数字(作回归对照, 偏离就说明改坏了):

| 脚本 | 期望 | 关键数字 |
|---|---|---|
| `run_validation.py` | 8/8 | SMIB 数值 CCT vs 等面积闭式差 **1e-6**; 潮流 9e-15; SEP 2.6e-11 |
| `run_matlab_xval.py` | 8/8 | Yred 1e-10; CUEP 1e-11; LEA 0.2274 vs 0.2275; REA 0.2433 vs 0.2434 |
| `test_config.py` | 5/5 | 9母线 SEP 2.6e-11; 39母线 SEP 1.7e-15 |

---

## 2. 日常交互: `cli.py` 的 5 个子命令

```bash
python cli.py show          # 打印当前配置摘要(不跑实验, 先看清楚要跑什么)
python cli.py list          # 列出全部可选 mode / case
python cli.py run           # 用 config.yaml 跑一次实验(含运行后自检)
python cli.py validate      # = P0 物理不变量 + 金标准(等价 run_validation 的验证段)
python cli.py xval          # = T3 与 MATLAB 逐层交叉验证
```

`list` 的真实输出:

```
可选 mode:
  - reduced_cct          reduced_numerical    reduced_region
  - spm_cct              spm_numerical        spm_region
  - two_machine_region_3d                     two_machine_region_3d_gfl
可选 case: case9_v2, case39_modified
```

---

## 3. 跑实验: 两个标准动作

### 动作 A — 9 母线, 全功能(能量法 LEA + 时域 REA)

```bash
python cli.py run
```

真实输出(节选, 去掉了配置摘要段):

```
---------- 运行 mode: reduced_cct ----------
[reduced_cct] LEA CCT(能量法) = 0.2274s   REA CCT(时域/有界判据) = 0.24s
[reduced_cct] 退出点 index = 199   CUEP 来源: closest-UEP (type-1, V>V(SEP))
[reduced_cct] 已保存图: ...\python_bcu\figures\reduced_cct_phaseplane.png
---------- 运行后自检 ----------
  postfault SEP 残差 : 1.04e-16  (PASS)
  能量法 CUEP        : type-1, V(CUEP)=2.1112, 候选 type-1 UEP=2
  LEA / REA CCT      : 0.2274 / 0.2433 s
  能量法保守 LEA<=REA: PASS
```

**你要记录的结果**: `LEA = 0.2274 s`, `REA = 0.2433 s`(mode 段与自检段现已一致, 均可信)。

### 动作 B — 39 母线, 能量法 LEA + 时域 REA

```bash
python cli.py run --case case39_modified --auto-params
```

真实输出(节选):

```
[config] 39 母线: faultline 由默认 [9,6] 自动切换为 [16,17](可用 --fault-line 覆盖).
---------- 运行 mode: reduced_cct ----------
[reduced_cct] LEA CCT(能量法) = 0.0847s   REA CCT(时域/有界判据) = 0.12s
[reduced_cct] 退出点 index = 199   CUEP 来源: closest-UEP (type-1, V>V(SEP))
---------- 运行后自检 ----------
  postfault SEP 残差 : 7.43e-09  (PASS)
  能量法 CUEP        : type-1, V(CUEP)=2.3732, 候选 type-1 UEP=8
  LEA / REA CCT      : 0.0847 / 0.1224 s  (ngen=10)
  能量法保守 LEA<=REA: PASS
```

**你要记录的结果**: 39 母线 `LEA = 0.0847 s`(能量法) / `REA = 0.1224 s`(时域精确)。
> ✅ 更新: 能量法 CUEP 已通用化到任意机数(`bcu_v2/cuep.py` closest-UEP), **39 母线也能出 LEA**(3 机对
> MATLAB 1e-11); 2026-09-01 该法回灌 v1 后, **mode 段与自检段一致**(退出点 199, 不再是 index=0 噪声)。

---

## 4. 怎么看懂输出

`cli.py run` 打印两段结果: **mode 段** + **运行后自检段**。

> ✅ 2026-09-01 重大更新: v1 的 6 个已知隐患已**就地修复**(见 `交接文档` 第 7 节), **mode 段现在也
> 正确了, 与自检段一致**。以前的两个"噪声陷阱"(9 母线 mode 段 LEA=0.057 错值、39 母线退出点=0)**已消失**。

| 段落 | 内容 | 现状 |
|---|---|---|
| `---------- 运行 mode ----------`(`[reduced_cct]` 行) | v1 实验模式过程输出(LEA/REA/退出点/CUEP来源/图) | ✅ 已修正, 可信 |
| `---------- 运行后自检 ----------` | SEP 残差 + 能量法 LEA/REA + 保守性检查(结构化汇总) | ✅ 可信 |

现在两段的 CCT 一致, 都可信:
- **9 母线**: LEA=0.2274 / REA≈0.24 s, 退出点 index=199, CUEP 来源=closest-UEP。
- **39 母线**: LEA=0.0847 / REA=0.1224 s, 退出点 index=199, CUEP 来源=closest-UEP。

> 历史注记: 修复前 mode 段用 v1 的 `find_mgp`(求错)得 9 母线 LEA=0.057、退出点=0, 那时需"只看自检段";
> 2026-09-01 v1 回灌修复(`find_mgp`→`closest-UEP`、退出点排除伪过零)后, mode 段与自检段一致, 无需再区分。

---

## 5. 现在能仿真什么 / 不能仿真什么(机数 × 模型矩阵)

从 `cli.py` 交互入口能实际跑到的:

| mode | 模型本质 | 发电机模型 | 9 母线(3机) | 39 母线(10机) |
|---|---|---|---|---|
| `reduced_cct` / `reduced_numerical` | 网络约简 | 经典二阶 | ✅ LEA + REA | ✅ **LEA + REA**(v2 通用 CUEP) |
| `reduced_region` | 稳定域(2D) | 经典二阶 | ✅ | ❌ 3 机专属 |
| `spm_cct` / `spm_numerical` | 结构保持模型 | 经典二阶 | ✅ | ❌ 未适配 |
| `spm_region` | 稳定域(SPM) | 经典二阶 | ✅ | ❌ 3 机专属 |
| `two_machine_region_3d(_gfl)` | 两机 3D 稳定域 | 经典二阶 | — | — (固定 2 机, 与电网无关) |

**要点:**
- 交互层面能跑的**发电机模型只有经典二阶**(`E'∠δ` 恒定)。这一点 v2 和 v1 相同。
- v2 相对 v1 的真实增量: **系统规模 9→39 母线** + **时域 CCT 精度**(事件驱动) + **可信度**(P0/T3)
  + **能量法 LEA 通用化到任意机数**(2026-08-31, 39 母线也能出 LEA)。
- **one-axis 同步机 / GFM 下垂 / ZIP 负荷 / 严格 DAE-SPM** 这些 v2 新写的模型**都还没接入
  `cli.py`**, 只能用 `python test_p2.py` / `test_spm_dae.py` 触发验证, 且都是**单机无穷大母线级**
  参考实现, 跑不了完整暂稳实验。

> 提醒: 39 母线的 H/Xd' 目前是**示例值**(`systems.NE39_H/NE39_XD1`), 定量 CCT 前务必用权威
> 参考(如 Pai 1989 附录)核对替换; 装配机理与 SEP/潮流不依赖其精确性。

---

## 6. 调参: 两种方式

### 方式一: 改 `config.yaml`(持久, 推荐用于稳定的实验设定)

打开 `config.yaml`, 只改 `key: value` 右边的值。核心参数:

| 参数 | 含义 |
|---|---|
| `mode` | 8 个实验模式之一(见第 5 节) |
| `case` | `case9_v2`(3机9母线) 或 `case39_modified`(10机39母线) |
| `auto_params` | `true` = 用案例内置动态参数(切 39 母线时设 true); `false` = 用下面手填的数组 |
| `m / damping_ratio / Pm / xd1 / E` | 每机参数, 长度必须 == 发电机台数(9母线=3, 否则校验报错) |
| `faultline / faultposition` | 故障支路 `[FromBus, ToBus]` 与位置 |
| `Tfault / Tunit / postfault_time` | 故障时长 / 积分步长 / 故障后仿真时长 |
| `region_grid` | 稳定域网格密度(每维点数) |

### 方式二: 命令行临时覆盖(一次性, 不动文件)

```bash
python cli.py run --mode reduced_region --grid 15
python cli.py run --case case39_modified --auto-params
python cli.py run --fault-line 8,9 --tunit 5e-4 --tfault 0.5
python cli.py run --set Pm=[0.9,1.3,0.95]        # 通用 key=value 覆盖(值走 JSON 解析)
```

常用覆盖参数: `--mode --case --auto-params --fault-line --tfault --tunit --grid --set`。
`--set` 可多次, 覆盖任意 config 键。

> 每次 `run` 会把配置快照存到 `results/snapshot_<时间戳>.json`(可追溯); 关掉设 `save_snapshot: false`。

---

## 7. 各阶段独立验证脚本(回归用)

改完代码要确认没破坏已验证行为时跑:

```bash
python run_validation.py   # P0 物理不变量 + 金标准(8/8)
python run_matlab_xval.py  # T3 与 MATLAB 交叉验证(8/8)
python test_fixes.py       # v1 五隐患修复(5/5)
python test_spm_dae.py     # 严格 DAE-SPM(3/3, DAE 7/7 vs v1 3/7)
python test_solvers.py     # scipy 求解器层(3/3)
python test_p2.py          # 建模扩展 39母线/ZIP/one-axis/GFM(6/6)
python test_config.py      # 配置系统 + T3 回归(5/5)
python test_cuep.py        # 通用能量法 CUEP + LEA(6/6): 9母线对 MATLAB 1e-11, 39母线出 LEA
```

---

## 8. 改代码前的规矩(硬约束)

1. **v1 已解冻(2026-09-01)**: 6 隐患已回灌修复、入口已移植。日常仍**优先在 `bcu_v2/` 改**; 确需改 v1
   须先跟负责人确认, 且**改 v1 后必须重跑 v2 全套**(v2 `import bcu_3m9b`, 会随 v1 变化)。
2. **MATLAB 平台只读**: 不写入 `../matlab_platform/`; T3 只加载其导出的 `verify/baseline_reduced.mat`。
3. **验证把关**: 每次改完, 先跑 `python run_validation.py && python run_matlab_xval.py &&
   python test_config.py`(改了 v1 再加 `python test_fixes.py && python test_cuep.py`), 确认"改进不破坏
   已验证行为", 再提交。
4. **注释/输出规范**: 注释 = 中文操作说明 + 半角标点; v1/v2 计算模块的程序输出 = 英文;
   cli/config 的中文 UI 可用中文。可运行脚本开头用 `SetConsoleOutputCP(65001)` 修 Windows 乱码。

---

## 9. 一分钟决策树(我现在该敲哪条命令?)

- 想确认平台没坏 → 第 1 节三条验证命令。
- 想看当前要跑什么 → `python cli.py show`。
- 想出 **9 母线** 的 CCT(能量法 + 时域) → `python cli.py run`, 读**自检段**的 LEA/REA。
- 想出 **39 母线** 的 CCT(能量法 + 时域) → `python cli.py run --case case39_modified --auto-params`,
  读自检段的 LEA / REA(能量法 LEA 已通用化, 约 4s)。
- 想扫参数出对比表 → 循环调 `python cli.py run --set <参数>=<值>`, 收集自检段数字。
- 想画稳定域(仅 3 机) → `python cli.py run --mode reduced_region`(图存到 `../python_bcu/figures/`)。

---

*本手册聚焦"当下怎么操作"。项目全貌/开发历程/后续方向见 `交接文档_HANDOVER_CN.md`;*
*模块级说明见 `README_CN.md`; 路线图见 `改进与提升路线图_CN.md`。*
