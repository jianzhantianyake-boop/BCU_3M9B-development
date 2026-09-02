# BCU 平台架构与使用文档 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在统一开发仓库中建立一份可追溯的中文项目文档，解释 Python/MATLAB 技术架构、已确认的 MATLAB 坐标问题、跨平台一致性测试证据和 v2/MATLAB 的基础操作方法。

**Architecture:** 新文档以平台分层、数据流、验证层和操作层组织，不复制八路径矩阵的全部内容，而是链接到现有权威矩阵和审计记录。所有技术结论绑定当前仓库提交、参考文件、命令、误差和限制；历史坐标混用结果与当前同坐标参考明确分开。

**Tech Stack:** Python 3.12、NumPy、SciPy、可选 Matplotlib、MATLAB R2024a、MATPOWER 7.1、PowerShell、Git。

**Spec:** 用户请求：详细总结当前 Python/MATLAB 技术架构和技术路径，明确 MATLAB bug 及证据，说明一致性测试，并提供 v2、MATLAB 和测试的基础使用说明。

## Global Constraints

- 只修改集成仓库，不修改三个来源仓库。
- 不把与 MATLAB 一致单独写成物理正确性证明。
- 不把历史 `3.3757/0.2053` 当作当前自足 SPM 验收目标。
- 明确区分已确认的 `deltacoi` 坐标混用与尚未证明执行的 `voltage_MGP` 疑似复制错误。
- 所有命令使用明确解释器和工作目录；失败、阻塞和不可比较状态不得用旧数字替代。
- 阶段 7 阻尼研究不执行，不在文档中生成阻尼实验结论。

---

### Task 1: 盘点现有文档和证据入口

**Files:**
- Read: `docs/01_项目现状与边界.md`
- Read: `docs/02_开发交接.md`
- Read: `docs/03_验证覆盖矩阵.md`
- Read: `docs/04_操作与回归指南.md`
- Read: `docs/05_开发路线图.md`
- Read: `docs/06_实验与结果规范.md`
- Read: `docs/provenance/2026-09-02_八路径总门禁.md`

- [x] **Step 1: 固定文档结构**

新文档必须包含：平台定位、目录和模块地图、Python v2 技术路径、MATLAB 技术路径、
跨平台数据契约、MATLAB bug 证据、测试分层、基础使用、常见失败和当前限制。

### Task 2: 编写权威架构与证据文档

**Files:**
- Create: `docs/07_平台技术架构与验证使用指南.md`
- Modify: `README.md`

**Interfaces:**
- References: `python_bcu_v2/bcu_v2/spm_energy.py`, `python_bcu_v2/bcu_v2/spm_dae.py`,
  `python_bcu_v2/bcu_v2/spm_cuep.py`, `python_bcu_v2/run_full_xval.py`。
- References: `matlab_platform/B3_MM/Cal_MM_Static_SPM.m`、`Cal_MM_CCT_SPM.m`、
  `Fun_Cal_Exitpoint_SPM.m`、`Fun_Cal_MGP_SPM.m`、`Fun_SEPfslove_SPM.m` 和
  `Fun_Cal_PotentialEnergy_SPM.m`。
- References: `validation/references/spm_cct_v1.json`、`spm_cct_v2.json` 和
  `docs/provenance/2026-09-02_八路径总门禁.md`。

- [x] **Step 1: 写平台架构章节**

说明 v1 兼容基线、v2 主开发平台、MATLAB 参考平台、MATPOWER 外部依赖、`.mat`/results
排除规则，以及 v2 通过 `../python_bcu` 复用 v1 包的边界。

- [x] **Step 2: 写技术路径章节**

分别给出 reduced、SPM、two-machine 的数据流；说明 SPM 的 fault1 DAE、MGP 连续追踪、
联合发电机—网络平衡、type-1 分类、五项势能和 LEA CCT 链路。

- [x] **Step 3: 写 MATLAB bug 证据章节**

记录 `postfault.CUEP_net_theta` 使用未限定 `deltacoi` 的代码表象、projected/raw 角度差、
正确物理网络残差、同坐标候选的 MATLAB 原生残差和能量结果。明确 `voltage_MGP =
theta_iter(no_MGP,:)` 目前只列为疑似问题，因为没有分支执行证据证明它影响了成功结果。

- [x] **Step 4: 写一致性测试章节**

按公式/不变量、Python 内部、MATLAB 原生、逐检查点、八路径总门禁和变异测试分层，列出
实际命令、检查数量、最大误差和不覆盖的范围。

- [x] **Step 5: 写基础使用章节**

给出 Python v2 的环境检查、单元测试、SPM 入口和统一报告命令；给出 MATLAB 的路径初始化、
静态初始化、CCT/SPM 顺序和 verify 诊断器命令；说明每一步的预期输出和 BLOCKED/FAILED 判定。

- [x] **Step 6: 更新 README 入口**

在 README 文档入口列表加入新指南和八路径总门禁记录链接，不在 README 重复完整矩阵。

### Task 3: 质量检查与版本提交

**Files:**
- Modify: `docs/修改日志.md`

- [x] **Step 1: 检查内部链接和过时断言**

运行：

```powershell
.\scripts\Test-Docs.ps1 -RepoRoot (Get-Location).Path
```

预期：`DOCS_OK`，且新文档中的路径全部可解析。

- [x] **Step 2: 检查 Markdown 结构**

运行：

```powershell
git diff --check
```

预期：无空白错误；使用 `rg` 搜索新文档中的 `3.3757`、`0.2053`、`7.57668`，确认每次出现
都带有“历史/坐标混用/不可作为当前目标”的限定。

- [x] **Step 3: 提交并验证**

```powershell
git add README.md docs/07_平台技术架构与验证使用指南.md docs/修改日志.md docs/superpowers/plans/2026-09-02-platform-architecture-and-usage-documentation.md
git commit -m "docs: document platform architecture and validation workflow"
git push -u origin main
git status --short --branch
git rev-list --left-right --count HEAD...origin/main
```

预期：工作区干净，`HEAD...origin/main` 为 `0 0`，作者仍为用户现有 Git 身份。
