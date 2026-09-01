% =========================================================================
% export_numerical.m -- 导出 reduced_numerical 的故障段末端切除态(代表性对比)
% numerical 完整为 pre(20s)+fault(0.24s)+post(80s)@1e-4(百万步); 其静态/CCT 已由
% reduced_cct(T3)覆盖, 此处取故障段末端(从 SEP 积分故障网络 0.24s)作代表性轨迹对比.
% =========================================================================
thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisDir);
addpath(projectRoot); setup_bcu_paths(projectRoot); cd(projectRoot);
load(fullfile(thisDir, 'baseline_reduced.mat'), 'preset', 'prefault', 'fault', 'Basevalue');

omegab = Basevalue.omegab;
ngen = numel(preset.m);
delta0 = prefault.SEP_delta(:);
omega0 = prefault.SEP_omegapu * omegab * ones(ngen, 1);

% 故障段: 从预故障 SEP 积分故障网络 0.24 s (numerical 的 Trecover-Tfault).
[theta, omega, thetac, ~, ~, cyc] = Fun_TrajIter_SRF(0.24, 1e-4, fault.Yred, preset, delta0, omega0, omegab);
theta_end  = theta(cyc, :);
omega_end  = omega(cyc, :);
thetac_end = thetac(cyc, :);

save(fullfile(thisDir, 'baseline_numerical.mat'), 'theta_end', 'omega_end', 'thetac_end', '-v7');
fprintf('EXPORT_NUM_OK theta_end=[%s]\n', num2str(theta_end, '%.6f '));
