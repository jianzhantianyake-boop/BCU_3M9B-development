% =========================================================================
% run_bcu_sweep.m —— 批量扫参驱动（一次扫一个参数，自动汇总 CCT）
%
% 【怎么用】只改下面「扫描定义」三行，然后命令窗口运行：run_bcu_sweep
%   sweep.param  : 要扫描的 bcu_config 字段名（字符串）
%   sweep.values : 元胞数组，每个元素是该字段的一组完整取值（可标量或向量）
%   sweep.fixed  : 每组都附加的固定覆盖（如缩短 Tfault 加速）；无则 struct()
%
% 【实现要点】Cal_MM_CCT 会 clear 工作区，普通 for 循环的累积变量会被清空。
%   本脚本把扫描状态存进 appdata(0)、参数经 bcu_override(persistent) 注入，
%   两者都不受 clear 影响，因此循环得以跨多次 clear 稳定推进，全程不改文件。
%
% 输出：results/sweep_<时间戳>.{csv,mat} + 一张 CCT-参数曲线图。
% =========================================================================

% ============================ 扫描定义 ==================================
sweep.param  = 'damping_ratio';
sweep.values = {[0.05;0.05;0.05], [0.1;0.1;0.1], [0.2;0.2;0.2], [0.3;0.3;0.3]};
sweep.fixed  = struct('Tfault', 0.8);          % 覆盖窗口够到 CCT(~0.25s)之后，加速扫描
sweep.label  = 'damping ratio d/m (uniform)';  % 图轴/标题用
% =======================================================================

thisFile = mfilename('fullpath');
projectRoot = fileparts(thisFile);
cd(projectRoot);
setup_bcu_paths(projectRoot);

N = numel(sweep.values);
setappdata(0, 'BCU_SW', sweep);                % 扫描定义（clear 不清）
setappdata(0, 'BCU_SW_RES', cell(N,1));        % 结果累积（clear 不清）
fprintf('==== 批量扫参开始：param=%s, 共 %d 组 ====\n', sweep.param, N);

% 注意：run(Cal_MM_CCT) 的 clear 会清掉循环变量 i/sw/ov/N。因此循环控制状态
% 一律存进 appdata(0)，并在 run 之后重新取回，绝不跨 run 直接引用工作区变量。
for i = 1:N
    setappdata(0, 'BCU_SW_I', i);              % 存当前索引
    sw = getappdata(0, 'BCU_SW');
    ov = sw.fixed;
    ov.(sw.param) = sw.values{i};
    bcu_override(ov);                          % 注入本组参数（persistent，clear 不清）
    fprintf('\n---- [%d/%d] %s = %s ----\n', i, numel(sw.values), sw.param, mat2str(sw.values{i}(:).'));
    run(fullfile('B3_MM','Cal_MM_CCT.m'));     % <== clear 在此；appdata/persistent 存活

    ii = getappdata(0, 'BCU_SW_I');            % run 后重取索引（i 已被 clear）
    sw = getappdata(0, 'BCU_SW');              % run 后重取扫描定义
    r  = getappdata(0, 'BCU_SW_RES');
    r{ii} = struct('idx', ii, 'value', sw.values{ii}, ...
                   'LEA', Critical.LEA.CCT, 'REA', Critical.REA.CCT, ...
                   'CUEP', postfault.CUEP_delta);
    setappdata(0, 'BCU_SW_RES', r);
    fprintf('[[SWEEP]] i=%d LEA=%.4f REA=%.4f\n', ii, r{ii}.LEA, r{ii}.REA);
end
bcu_override('clear');                         % 扫参结束务必清除覆盖，避免影响后续单次运行

% ============================ 汇总与输出 ================================
r  = getappdata(0, 'BCU_SW_RES');
sw = getappdata(0, 'BCU_SW');
N  = numel(sw.values);                         % 重建 N（循环里已被 clear）
LEA = zeros(N,1); REA = zeros(N,1); xlab = cell(N,1);
for i = 1:N
    LEA(i) = r{i}.LEA; REA(i) = r{i}.REA;
    xlab{i} = mat2str(r{i}.value(:).');
end

projectRoot = fileparts(mfilename('fullpath'));   % 重建（循环里已被 clear）
resdir = fullfile(projectRoot, 'results');
if ~isfolder(resdir); mkdir(resdir); end
stamp = datestr(now, 'yyyymmdd_HHMMSS');
T = table((1:N)', xlab, LEA, REA, 'VariableNames', {'idx','value','LEA_CCT','REA_CCT'});
writetable(T, fullfile(resdir, ['sweep_' stamp '.csv']));
save(fullfile(resdir, ['sweep_' stamp '.mat']), 'sw', 'r', 'LEA', 'REA');

fprintf('\n================ 扫参结果：%s ================\n', sw.label);
for i = 1:N
    fprintf('  [%d] %-20s LEA-CCT=%.4f  REA-CCT=%.4f s\n', i, xlab{i}, LEA(i), REA(i));
end
fprintf('  结果已存 results/sweep_%s.{csv,mat}\n', stamp);

figure('Name','CCT sweep');
plot(1:N, LEA, '-o', 1:N, REA, '-s', 'LineWidth', 1.6); grid on; hold on;
set(gca, 'XTick', 1:N, 'XTickLabel', xlab);
try; xtickangle(25); catch; end
xlabel(sw.label, 'Interpreter', 'none'); ylabel('CCT (s)');
legend('LEA-CCT (energy)', 'REA-CCT (real)', 'Location', 'best');
title('CCT vs sweep parameter');
