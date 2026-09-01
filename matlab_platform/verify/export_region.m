% =========================================================================
% export_region.m -- 导出 reduced_region 路径的 MATLAB 参考(平衡点集合)
% 用法: matlab -batch "run('matlab_platform/verify/export_region.m')"
% 不改 B3_MM; 运行原 mode(抑制画图), 抓取 base workspace 的 ep_set, 精简存 .mat.
% =========================================================================
thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisDir);
addpath(projectRoot); setup_bcu_paths(projectRoot);
set(0, 'DefaultFigureVisible', 'off');   % 无界面: 抑制图窗
cd(projectRoot);

% reduced_region 链路: 先静态初始化, 再稳定域 EP 搜索(同 run_bcu.m).
run(fullfile('B3_MM', 'Cal_MM_Static.m'));
run(fullfile('B3_MM', 'Statable_Region.m'));

% 提取平衡点集合: xep=[d2c,d3c], flag=非负特征值个数(0=SEP,1=type-1 UEP,...).
nep = numel(ep_set);
xeps = zeros(nep, 2);
flags = zeros(nep, 1);
for i = 1:nep
    xeps(i, :) = ep_set(i).xep(:).';
    flags(i) = ep_set(i).flag;
end

outfile = fullfile(thisDir, 'baseline_region.mat');
save(outfile, 'xeps', 'flags', '-v7');
fprintf('EXPORT_REGION_OK nep=%d file=%s\n', nep, outfile);
