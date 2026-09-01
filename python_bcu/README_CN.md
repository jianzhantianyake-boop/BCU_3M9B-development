# Python BCU 实验仿真平台

这是 `BCU_3M9B-main` 的独立 Python 实验层。目标是让学习者可以在不启动
MATLAB 的情况下，逐步观察潮流、导纳矩阵、网络约简、SEP、MGP、CUEP、能量
和暂态轨迹。

当前实现依赖 NumPy 与 SciPy（SciPy 用于稳健求解器：SPM 负荷代数解、CUEP/closest-UEP
的平衡点求根）；Matplotlib 只由 `plotting.py` 按需导入。

> 2026-09-01 更新：v1 已解除冻结，6 个已知隐患全部就地修复（solve_sep 副作用/CUEP 污染/
> 退出点伪过零/轻阻尼判据/SPM 代数解冷启动/find_mgp 求错），并新增**配置驱动入口** `cli.py`
> + `config.yaml`（详见下方"入口"）。能量法 LEA CCT 已由 closest-UEP（`bcu_3m9b/cuep.py`）
> 修正，与 MATLAB 平台一致到 ~1e-11。

## 入口（两种）

- **交互菜单**（探索用）：`python main.py`（8 模式菜单 + 工具项）。
- **配置驱动**（推荐，可复现）：编辑 `config.yaml` 后 `python cli.py run`（含运行后自检）；
  `python cli.py show / list / validate`；命令行临时覆盖，如
  `python cli.py run --mode reduced_region --grid 15` 或 `--set Pm=[0.9,1.3,0.95]`。

## 最小运行

```powershell
cd "C:\Users\WangS\Desktop\26年第二学期学习文件\论文文献\literature_ntu_chatgpt_only\代码\BCU_3M9B-main\python_bcu"
& "C:\Users\WangS\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" scripts\run_static.py
& "C:\Users\WangS\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" scripts\run_experiment.py
& "C:\Users\WangS\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" tests\smoke_test.py
```

如果电脑已安装 Python 和 NumPy，可把命令中的绝对 Python 路径替换成
`python`。不建议把项目路径写入环境变量，脚本本身会把 `python_bcu` 加入
模块搜索路径。

## 推荐运行顺序

1. `scripts/run_two_machine.py`：看两机摆动方程、平衡点和特征值分类。
2. `scripts/run_static.py`：查看 9-bus 潮流、发电机顺序、Yred 和 SEP。
3. `scripts/run_spm.py`：查看结构保持模型的负荷母线代数方程。
4. `scripts/run_experiment.py`：运行故障轨迹、退出点、MGP 候选、能量 CCT 和
   时域网格 CCT。

## 平台对象

`CaseData` 保留 MATPOWER 的 bus/gen/branch 原始列顺序；`PFData` 对应原工程
`Fun_ResultBack` 的项目结构；`Preset` 保存 m、d、Pm、E 和故障设置；
`NetworkState` 保存 Yfull、Yred、分块矩阵和 SEP；`Trajectory` 保存固定步长
轨迹；`StaticResult` 汇总 `Cal_MM_Static` 的全部主要输出。

角度统一使用 rad（仅打印潮流母线电压时显示 deg），角速度使用 rad/s，功率
和导纳使用 pu，潮流输入输出的 P/Q 使用 MW/MVAr。发电机顺序默认是母线
1、2、3；故障线路默认是 9--6，故障母线默认是 9。

## 研究实现与平台实现的边界

Python 版本保持原项目的主要公式和数据顺序，但把 MATLAB 的手写迭代、
`fsolve`、`ode15s` 和交互式 `input` 改成了可批处理的 NumPy/SciPy 接口。SPM 使用
“每个显式时间步先解负荷代数方程（scipy 稳健求解），再推进发电机”的可运行近似；
CUEP 改用 closest-UEP（`bcu_3m9b/cuep.py`，`find_mgp` 降为诊断/回退）。**核心链路
（潮流/Yred/SEP/CUEP/CCT）已通过 v2 的 T3 与 MATLAB 平台逐层交叉验证到 ~1e-10~1e-11**
（v2 `import bcu_3m9b`，验证的就是本层代码）。绘图、稳定域网格和大规模参数扫描可以
在此核心层之上继续扩展。

MATPOWER 官方源码未复制、未修改。Python 案例数据仅转写了项目自有的
`case9_v2.m` 和 `case39_modified.m`。
