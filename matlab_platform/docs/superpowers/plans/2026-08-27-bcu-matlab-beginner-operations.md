# BCU MATLAB 新手操作入口与教学注释 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a safe single-experiment MATLAB launcher, a beginner operating guide, and detailed operation-style Chinese comments without changing any BCU or MATPOWER algorithm.

**Architecture:** A root-level base-workspace script selects one existing entrypoint and delegates to it. Documentation describes dependencies and parameter edits rather than introducing configuration overrides. Existing project-owned MATLAB files receive comments only; their executable statements remain byte-for-byte unchanged.

**Tech Stack:** MATLAB scripts/functions, bundled MATPOWER 7.1, Markdown documentation, PowerShell read-only/static validation.

**Spec:** `docs/superpowers/specs/2026-08-27-bcu-matlab-beginner-operations-design.md`

## Global Constraints

- The controller is a MATLAB script executed in the base workspace; do not convert it to a function.
- Do not change equations, parameter defaults, function signatures, solver calls, algorithmic criteria, or MATPOWER official source.
- The controller supports only one manually selected experiment at a time; no sweeps, exports, retries, or hidden parameter overrides.
- Every changed project-owned MATLAB header contains 使用方法、参数、返回或工作区结果、步骤、单位、前置条件 and 研究/验证边界.
- MATLAB native execution is `BLOCKED` by the known startup failure; report static verification separately from E2E.

---

### Task 1: Add the base-workspace experiment launcher

**Files:**

- Create: `run_bcu_beginner.m`
- Read: `setup_bcu_paths.m`, `B3_MM/Cal_MM_CCT.m`, `B3_MM/Cal_MM_CCT_SPM.m`, `B3_MM/NumSim_MM_Gridframe.m`, `B3_MM/NumSim_MM_Gridframe_SPM.m`

**Interfaces:**

- Consumes: one user-edited string variable `EXPERIMENT_MODE` in the script.
- Produces: the original scripts' base-workspace structures, figures and command-window messages; no new result data structure.

- [ ] **Step 1: Confirm the required entrypoints exist**

Run a read-only path check for `setup_bcu_paths.m`, the four CCT/numerical entry scripts and the four stable-region scripts. Stop if any named entrypoint is missing.

- [ ] **Step 2: Create a commented selector script**

Implement the modes `reduced_cct`, `reduced_numerical`, `spm_cct`, `spm_numerical`, `reduced_region`, `spm_region`, `two_machine_region_3d`, and `two_machine_region_3d_gfl`. Use `run(fullfile('B3_MM', ...))` after switching MATLAB's current folder to the controller's folder, so later calls remain valid after CCT scripts issue `clear`.

- [ ] **Step 3: Make numerical modes prepare CUEP data**

For both numerical modes, execute the matching CCT script before the matching `NumSim` script. Do not add synthetic `CUEP_delta` values or alter numerical scripts.

- [ ] **Step 4: Keep unsupported plotting helpers out of automatic dispatch**

Explain in comments that `Plot_3Dstate.m` needs its named MAT file and `plottmp.m` needs `Group`; do not add automatic modes that hide those prerequisites.

### Task 2: Write the beginner operating guide

**Files:**

- Create: `docs/BCU_MATLAB新手操作指南.md`
- Read: `README.md`, `docs/BCU学习与运行指南.md`, `run_bcu_beginner.m`

**Interfaces:**

- Consumes: the project folder, MATLAB, bundled MATPOWER and the controller mode names.
- Produces: reproducible manual running instructions and parameter-change locations.

- [ ] **Step 1: Explain first startup and path audit**

Document `paths = setup_bcu_paths();`, `which runpf -all`, `which case9_v2 -all`, and `which Fun_ResultBack -all`. State that all results are blocked until MATLAB starts without the known settings-plugin failure.

- [ ] **Step 2: Document every single-experiment mode**

For each launcher mode, state prerequisite, exact original scripts executed, main workspace outputs, automatic figures, and what a newcomer should inspect before interpreting results.

- [ ] **Step 3: Add parameter-edit map**

Map generator/base parameters, `PathEnergyCal`, case selection, fault line/position, CCT fault duration, numerical time settings and solver-related settings to their existing source locations. State whether the file must be re-run from a clean workspace after an edit.

- [ ] **Step 4: Add a drawing and troubleshooting chapter**

Document automatic CCT/numerical figures, `Plot_3Dstate.m` data-file prerequisite, `plottmp.m` workspace prerequisite, `CUEP_delta` error recovery, base-workspace dependency, and warning against treating exploratory stable-region scripts as global stability proofs.

### Task 3: Convert project-owned MATLAB headers to operation-style Chinese comments

**Files:**

- Modify: `Y.m`, `setup_bcu_paths.m`
- Modify: all `B3_MM/*.m` files
- Modify: `C1_Matpower/matpower7.1/Fun_ResultBack.m`, `Fun_ResultSaved.m`, `ResultSaved.m`, `data/case9_v2.m`, `data/case39_modified.m`

**Interfaces:**

- Consumes: existing file purpose, function signature, script dependencies and base-workspace names.
- Produces: comments only; every executable statement and function signature remains unchanged.

- [ ] **Step 1: Annotate entry and initialization scripts**

Describe how to run `Y.m`, `setup_bcu_paths.m`, `Cal_MM_Static.m`, and `Cal_MM_Static_SPM.m`; identify `preset`, `prefault`, `fault`, `postfault`, `Basevalue`, paths and units.

- [ ] **Step 2: Annotate CCT and numerical simulation scripts**

Describe CCT preparation, CUEP dependency, three-stage time segments, time units, generated figures and the distinction between energy CCT and trajectory CCT.

- [ ] **Step 3: Annotate B3 helper functions by signature and model family**

For every helper, make the input/output relation visible from its current function signature; distinguish reduced-model and `_SPM` functions, state vectors, matrices and units.

- [ ] **Step 4: Annotate stable-region and plotting scripts**

State standalone versus workspace/data-file prerequisites, figure side effects, dimensions, and the exploratory nature of the result.

- [ ] **Step 5: Annotate project-owned MATPOWER bridge files and cases**

Document their call path, inputs/outputs and distinction from untouched official MATPOWER files.

### Task 4: Update the modification log and run static verification

**Files:**

- Modify: `docs/修改日志.md`
- Verify: all files listed in Tasks 1–3

**Interfaces:**

- Consumes: final file list and static read-only checks.
- Produces: an auditable modification record and a two-level validation statement.

- [ ] **Step 1: Append a dated change record**

List created files, comment-only files, algorithm-impact statement, static checks and the still-blocked MATLAB E2E state.

- [ ] **Step 2: Verify paths, comment-field coverage and excluded source boundary**

Check every controller target exists, every selected project-owned header contains the required operation fields, and no extra file under the official MATPOWER library changed.

- [ ] **Step 3: Review diff-level behavior boundary**

Confirm that all modified legacy MATLAB files differ only in `%` comment lines; if any executable line differs, stop and investigate before reporting completion.

