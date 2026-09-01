function paths = setup_bcu_paths(projectRoot)
% =========================================================================
% 教学操作说明：配置 BCU_3M9B 项目的可移植 MATLAB 搜索路径。
% 使用方法：
%   推荐在任何实验前执行：
%       paths = setup_bcu_paths();
%   如果项目被复制到其他位置，可显式指定根目录：
%       paths = setup_bcu_paths('C:\...\BCU_3M9B-main');
%   成功后再运行 run_bcu_beginner.m 或 B3_MM 中的原始实验入口。
% 参数：
%   projectRoot（可选，char 或 string）：项目根目录。省略或为空时，自动
%   取本文件所在目录；该目录必须包含 B3_MM 和 C1_Matpower/matpower7.1。
% 返回：
%   paths（结构体）：root、b3mm、matpowerRoot、matpowerLib、matpowerData、
%   mostLib、optModelLib、mipsLib、mpTestLib 等路径。函数同时调用 addpath。
% 步骤：
%   1. 确定项目根目录。2. 逐一检查项目自有目录是否存在。3. 将 B3_MM、
%   MATPOWER 主目录及所需库目录加入当前 MATLAB 会话的搜索路径。
% 单位：
%   所有输入/输出都是文件系统路径字符串，无物理单位。
% 前置条件：
%   项目目录结构必须完整；若 MATLAB 已经从其他版本 MATPOWER 加载同名函数，
%   请先检查 path 冲突，不要仅凭“未报错”判断加载的是本项目版本。
% 研究与验证边界：
%   本函数只处理路径，不改变 BCU 方程、网络矩阵、初值、数值容差或实验结果。
% =========================================================================

if nargin < 1 || isempty(projectRoot)
    thisFile = mfilename('fullpath');
    projectRoot = fileparts(thisFile);
end

projectRoot = char(projectRoot);
if ~isfolder(projectRoot)
    error('BCU:PathNotFound', '项目根目录不存在：%s', projectRoot);
end

paths.root = projectRoot;
paths.b3mm = fullfile(projectRoot, 'B3_MM');
paths.matpowerRoot = fullfile(projectRoot, 'C1_Matpower', 'matpower7.1');
paths.matpowerLib = fullfile(paths.matpowerRoot, 'lib');
paths.matpowerData = fullfile(paths.matpowerRoot, 'data');
paths.mostLib = fullfile(paths.matpowerRoot, 'most', 'lib');
paths.optModelLib = fullfile(paths.matpowerRoot, 'mp-opt-model', 'lib');
paths.mipsLib = fullfile(paths.matpowerRoot, 'mips', 'lib');
paths.mpTestLib = fullfile(paths.matpowerRoot, 'mptest', 'lib');

requiredDirectories = {paths.b3mm, paths.matpowerRoot, paths.matpowerLib, ...
    paths.matpowerData, paths.mostLib, paths.optModelLib, paths.mipsLib, ...
    paths.mpTestLib};
for k = 1:numel(requiredDirectories)
    directory = requiredDirectories{k};
    if ~isfolder(directory)
        error('BCU:PathNotFound', '项目依赖目录不存在：%s', directory);
    end
    addpath(directory);
end

% 项目自有的 MATPOWER 接口文件位于 matpower7.1 根目录。
addpath(paths.matpowerRoot);
% 项目根：确保 setup_bcu_paths / bcu_config / run_bcu 等根目录脚本在
% run() 切换当前目录（cd 到 B3_MM）后仍可被搜索到。
addpath(paths.root);

% 预热 MATPOWER 的 osqp 能力检测缓存。本机未编译 osqp_mex，而 Cal_MM_Static
% 会用 addpath(genpath(matpower)) 递归引入 @osqp 接口类；若等到那时才首次探测，
% have_feature 会在实例化 osqp 时报错（间歇性致命）。此处 osqp 接口尚不可见，
% 先探测一次让结果缓存为 0，之后不再触发。探测失败也无害（try/catch 吞掉）。
try
    have_feature('osqp');
catch
end
end
