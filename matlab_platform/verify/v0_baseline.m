% =========================================================================
% v0_baseline.m —— 网络约简模型端到端基准 + 第 2 层内部一致性自检
%
% 目的：
%   1. 确认 Cal_MM_Static + Cal_MM_CCT 在本机端到端跑通（冒烟测试）。
%   2. 抓取关键数值（SEP / CUEP / MGP / exit point / LEA-CCT / REA-CCT）
%      保存为 baseline_reduced.mat，作为后续交叉验证与回归的锚点。
%   3. 就地做第 2 层可算的一致性自检并以 [[VERIFY]] 前缀打印，便于抓取。
%
% 运行方式（二选一）：
%   - VSCode / MATLAB：把当前目录切到 matlab_platform，运行本文件。
%   - 命令行：matlab -batch "run('verify/v0_baseline.m')"（cwd=matlab_platform）
%
% 注意：Cal_MM_CCT 会 clear base workspace 并自动调用 Cal_MM_Static，
%       因此本脚本不在其之前保留任何依赖变量。
% =========================================================================

% 自定位：本文件在 <projectRoot>/verify/，父目录即项目根。
% 这样无论用命令行 run 还是 VSCode 的 Run 按钮（都会 cd 到本文件目录）都成立。
thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisDir);
addpath(projectRoot);       % 先让 setup_bcu_paths 可见
setup_bcu_paths(projectRoot);

tic;                        % 无参全局计时器：Cal_MM_CCT 的 clear 不会清除它
Cal_MM_CCT;                 % 端到端主流程（内部 clear + Cal_MM_Static）
elapsed = toc;              % 仅供参考；若内部函数调用过 tic 则此值偏小

% ---- 保存完整基准工作区 ----
% 注意：Cal_MM_CCT 已 clear，thisDir 变量已丢失；用 mfilename 现算存到 verify/。
outMat = fullfile(fileparts(mfilename('fullpath')), 'baseline_reduced.mat');
save(outMat);

% ============================ 第 2 层自检 ================================
fprintf('\n============ [[VERIFY]] v0 baseline 自检开始 ============\n');
fprintf('[[VERIFY]] 端到端运行耗时 = %.1f s\n', elapsed);

% --- 检查 1：SEP / CUEP 平衡残差应接近 0 ---
sep_pre  = norm(prefault.SEP_Perr);
sep_post = norm(postfault.SEP_Perr);
cuep_err = norm(postfault.CUEP_Perr);
fprintf('[[VERIFY]] prefault  SEP 残差 |Perr| = %.3e\n', sep_pre);
fprintf('[[VERIFY]] postfault SEP 残差 |Perr| = %.3e\n', sep_post);
fprintf('[[VERIFY]] postfault CUEP 残差 |Perr| = %.3e\n', cuep_err);
pass_equ = (sep_pre < 1e-6) && (sep_post < 1e-6) && (cuep_err < 1e-6);
fprintf('[[VERIFY]] 平衡残差检查: %s (阈值 1e-6)\n', tf2str(pass_equ));

% --- 检查 2：CUEP 必须显著区别于 SEP（否则找错了点）---
d_cuep_sep = norm(postfault.CUEP_delta - postfault.SEP_delta);
fprintf('[[VERIFY]] |CUEP - SEP| = %.4f (应 > 0.1)\n', d_cuep_sep);

% --- 检查 3：能量法 CCT 应保守，即 LEA-CCT <= REA-CCT ---
lea = Critical.LEA.CCT;
rea = Critical.REA.CCT;
fprintf('[[VERIFY]] LEA-CCT (能量直接法) = %.6f s\n', lea);
fprintf('[[VERIFY]] REA-CCT (逐步仿真真值) = %.6f s\n', rea);
fprintf('[[VERIFY]] CCT 相对差 = %.2f%%\n', 100*(rea-lea)/rea);
pass_cct = lea <= rea + 1e-9;
fprintf('[[VERIFY]] 能量法保守性 (LEA<=REA): %s\n', tf2str(pass_cct));

% --- 关键几何量打印（供与上游/ Python 交叉验证）---
fprintf('[[VERIFY]] prefault.SEP_delta  = %s\n', vec2str(prefault.SEP_delta));
fprintf('[[VERIFY]] postfault.SEP_delta = %s\n', vec2str(postfault.SEP_delta));
fprintf('[[VERIFY]] postfault.CUEP_delta= %s\n', vec2str(postfault.CUEP_delta));
fprintf('[[VERIFY]] MGP.thetac_MGP      = %s\n', vec2str(MGP.thetac_MGP));
fprintf('[[VERIFY]] escape.thetac       = %s\n', vec2str(escape.thetac));
fprintf('[[VERIFY]] baseline 已保存: %s\n', outMat);
fprintf('============ [[VERIFY]] v0 baseline 自检结束 ============\n');

% ------------------------- 局部辅助函数 --------------------------------
function s = tf2str(b)
    if b, s = 'PASS'; else, s = '**FAIL**'; end
end
function s = vec2str(v)
    s = ['[', sprintf('%.4f ', v(:).'), ']'];
end
