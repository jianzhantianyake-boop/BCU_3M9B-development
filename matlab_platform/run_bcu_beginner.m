% =========================================================================
% 单次实验控制脚本：BCU_3M9B MATLAB 平台新手入口
%
% 使用方法：
%   1. 在 MATLAB 中打开本文件所在的 BCU_3M9B-main 文件夹。
%   2. 只修改下方的 EXPERIMENT_MODE 一行，选择一次要运行的实验。
%   3. 点击“运行”。一次运行只执行一个实验链；本脚本不会批量扫参、
%      自动保存图片，也不会改写任何模型参数。
%   4. 想开始另一组实验时，先在命令窗口执行 clear; close all;，再改
%      EXPERIMENT_MODE 后重新运行本脚本。这样可避免 base workspace 中
%      的变量被不同模型混用。
%
% 参数：
%   EXPERIMENT_MODE  字符串。可选值及其运行顺序如下：
%   "reduced_cct"             网络约简模型：初始化 -> CCT / CUEP 计算。
%   "reduced_numerical"       网络约简模型：初始化 -> CCT / CUEP -> 数值轨迹。
%   "reduced_region"          网络约简模型：初始化 -> 二维稳定区域搜索。
%   "spm_cct"                 结构保持模型：初始化 -> CCT / CUEP 计算。
%   "spm_numerical"           结构保持模型：初始化 -> CCT / CUEP -> 数值轨迹。
%   "spm_region"              结构保持模型：初始化 -> 稳定区域搜索。
%   "two_machine_region_3d"   独立两机三维稳定区域示例。
%   "two_machine_region_3d_gfl" 独立的 GFL 相关三维探索示例。
%
% 返回 / 工作区结果：
%   本脚本本身不返回函数值。被调用的原始脚本会在 MATLAB base workspace
%   生成 preset、prefault、fault、postfault、CUEP 等变量，并按原实现
%   打开图窗。变量名和图窗数量随实验模式变化。
%
% 步骤：
%   1. 自动定位项目根目录并调用 setup_bcu_paths，避免依赖当前工作目录。
%   2. 根据 EXPERIMENT_MODE 调用原始初始化、CCT、轨迹或稳定区域脚本。
%   3. 不在这里覆盖原始脚本中的参数；需要改参数时，请按
%      docs/BCU_MATLAB新手操作指南.md 的“参数修改”部分逐项修改。
%
% 单位：
%   原模型通常采用角度 rad、角速度 rad/s、时间 s、功率/电压/导纳 pu。
%   本控制脚本不进行任何单位换算。
%
% 前置条件：
%   - MATLAB 能正常启动，且具备 Optimization Toolbox（fsolve）。
%   - 项目内的 ode78 与 MATPOWER 7.1 文件完整。
%   - CCT / 数值轨迹模式会执行耗时计算；不要在同一次运行中同时打开
%     多个模式，也不要在计算中手动清空工作区。
%
% 研究边界：
%   本文件仅提供安全的单次调度，不改变任何 BCU 方程、判据、参数或
%   原始绘图逻辑。当前环境尚未完成 MATLAB 原生端到端验证；请以实际
%   MATLAB 运行结果为准，并记录 MATLAB 版本、工具箱和参数设置。
% =========================================================================

% --- 新手只需要修改这一行；每次只能选择一个模式。---
EXPERIMENT_MODE = "two_machine_region_3d";

% --- 以下调度代码通常不要修改。---
thisFile = mfilename('fullpath');
projectRoot = fileparts(thisFile);
if ~isfolder(projectRoot)
    error('BCU:ProjectNotFound', '无法定位项目根目录：%s', projectRoot);
end

% 统一从项目根目录运行，避免原始脚本使用相对路径时找不到文件。
cd(projectRoot);
projectPaths = setup_bcu_paths(projectRoot); %#ok<NASGU>

fprintf('\nBCU 单次实验控制器已启动。模式：%s\n', EXPERIMENT_MODE);
fprintf('项目根目录：%s\n\n', projectRoot);

% 注意：Cal_MM_CCT*.m 原始脚本会 clear 工作区。因此“轨迹”模式中的
% 两个 run 必须保持如下顺序，且不要在它们之间插入依赖临时变量的代码。
switch EXPERIMENT_MODE
    case "reduced_cct"
        run(fullfile('B3_MM', 'Cal_MM_CCT.m'));

    case "reduced_numerical"
        run(fullfile('B3_MM', 'Cal_MM_CCT.m'));
        run(fullfile('B3_MM', 'NumSim_MM_Gridframe.m'));

    case "reduced_region"
        run(fullfile('B3_MM', 'Cal_MM_Static.m'));
        run(fullfile('B3_MM', 'Statable_Region.m'));

    case "spm_cct"
        run(fullfile('B3_MM', 'Cal_MM_CCT_SPM.m'));

    case "spm_numerical"
        run(fullfile('B3_MM', 'Cal_MM_CCT_SPM.m'));
        run(fullfile('B3_MM', 'NumSim_MM_Gridframe_SPM.m'));

    case "spm_region"
        run(fullfile('B3_MM', 'Cal_MM_Static_SPM.m'));
        run(fullfile('B3_MM', 'Statable_Region_SPM.m'));

    case "two_machine_region_3d"
        % 该脚本自带参数并会 clear；不依赖 9-bus 初始化结果。
        run(fullfile('B3_MM', 'Statable_Region_3D.m'));

    case "two_machine_region_3d_gfl"
        % 这是独立的探索性示例；“GFL”是原脚本命名，不代表已完成
        % 面向所有 GFL 场景的验证。
        run(fullfile('B3_MM', 'Statable_Region_3D_GFL.m'));

    otherwise
        error('BCU:UnknownExperimentMode', ...
            ['未知 EXPERIMENT_MODE：%s。请从 reduced_cct、' ...
             'reduced_numerical、reduced_region、spm_cct、' ...
             'spm_numerical、spm_region、two_machine_region_3d、' ...
             'two_machine_region_3d_gfl 中选择一个。'], EXPERIMENT_MODE);
end

fprintf('\n原始实验脚本已返回。请先检查命令窗口警告和图窗，再决定是否保存结果。\n');

