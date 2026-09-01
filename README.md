# BCU_3M9B development repository

这是 BCU_3M9B 平台的统一开发仓库：以 `python_bcu_v2/` 为主开发平台，以
`matlab_platform/` 为只读参考方程，以 `python_bcu/` 为兼容基线。三个目录是从
已记录的来源提交建立的干净快照；本仓库不会自动回写来源仓库。

## 文档入口

- [项目现状与边界](docs/01_项目现状与边界.md)
- [开发交接](docs/02_开发交接.md)
- [验证覆盖矩阵](docs/03_验证覆盖矩阵.md)
- [操作与回归指南](docs/04_操作与回归指南.md)
- [开发路线图](docs/05_开发路线图.md)
- [实验与结果规范](docs/06_实验与结果规范.md)
- [修改日志](docs/修改日志.md)
- [来源审计](docs/provenance/2026-09-01_接手来源清单.md)

## 维护原则

来源提交、工作区状态、文件哈希和排除理由记录在根目录
`SOURCE_MANIFEST.csv` 及 `docs/provenance/`。大型 MATPOWER 依赖、MATLAB `.mat`
工作区、运行结果和缓存不进入 Git，需要时按操作指南在本地引导或重生。

所有验证结果必须区分物理不变量验证、Python 内部验证、MATLAB 交叉验证、近似实现、
未验证和阻塞状态。“与 MATLAB 一致”本身不是物理正确性的证明；外部
`E_critical` 不能被冒充为 SPM 自足结果。

## 最短维护流程

1. 运行 `scripts/audit_environment.ps1` 记录解释器、MATLAB、MATPOWER 和 Git 状态。
2. 按 [操作与回归指南](docs/04_操作与回归指南.md) 运行统一验证入口。
3. 更新代码、验证报告和对应权威文档；失败或无法运行时保留 `BLOCKED`/`FAILED`。
4. 提交前检查 `SOURCE_MANIFEST.csv` 未被意外改写，且没有 `.mat`、结果或凭据进入暂存区。

本仓库的 GitHub 远程必须保持私有；远程创建、身份验证和推送遵循交接审计中的安全边界。

