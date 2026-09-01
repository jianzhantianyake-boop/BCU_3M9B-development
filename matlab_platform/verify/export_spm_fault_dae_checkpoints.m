% =========================================================================
% export_spm_fault_dae_checkpoints.m -- 导出与 fault1 DAE 一致的 SPM 检查点
%
% 旧的 fault.traj 字段来自 Cal_MM_CCT_SPM 的故障后网络校正，不能作为
% fault1 代数方程的参考。本脚本只读运行 Cal_MM_Static_SPM，并复用
% F_3M9B_SP_DAE + ode15s 生成真正 fault1 轨迹；不修改 B3_MM 文件。
% =========================================================================
thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisDir);
addpath(projectRoot); setup_bcu_paths(projectRoot);
set(0, 'DefaultFigureVisible', 'off');
cd(projectRoot);

% Cal_MM_Static_SPM 会 clear 当前 base workspace；路径变量在运行后重新解析。
run(fullfile('B3_MM', 'Cal_MM_Static_SPM.m'));
thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisDir);

checkpointTimes = [0, 0.1, 0.2, 0.2053, 0.21, 0.25, 0.3];
Tunit = 1e-4;
Tfault = max(checkpointTimes);
ngen = preset.ngen;
nbus = preset.nbus;
nnet = nbus - ngen;
faultbus = fault.faultbus;
placeholder = faultbus - ngen;

delta0 = prefault.SEP_delta(:);
omega0 = prefault.SEP_omegapu * Basevalue.omegab * ones(ngen, 1);
delta_net0 = prefault.net_delta(:);
voltage_net0 = prefault.net_voltage(:);
delta_net0(placeholder) = [];
voltage_net0(placeholder) = [];
% 保留原 NumSim_MM_Gridframe_SPM 的六槽状态布局：第六槽是删去故障母线的占位。
delta_net0(nnet) = 0;
voltage_net0(nnet) = 0;
[delta_net_s, V_net_s, flag_iter, ~, err] = ...
    Fun_AEiteration_SPM(delta_net0, voltage_net0, delta0, preset, Basevalue, ...
                        "fault1", 1e5, 1e-10);
if flag_iter ~= 1 || any(~isfinite([delta_net_s(:); V_net_s(:)]))
    error('SPM_FAULT_EXPORT_INIT: fault1 algebraic initialization failed (err=%g)', err);
end

M = diag([ones(ngen,1); ones(12,1)*1e-10; ones(ngen,1)]);
options = odeset('Mass', M, 'RelTol', 1e-10, ...
    'AbsTol', [1e-8*ones(1,ngen), 1e-12*ones(1,12), 1e-8*ones(1,ngen)]);
x0 = [delta0; delta_net_s(:); V_net_s(:); omega0];
tgrid = (0:Tunit:Tfault).';
[tout, xall] = ode15s(@(t,x)F_3M9B_SP_DAE(x, "fault1"), tgrid, x0, options);

idx = zeros(numel(checkpointTimes), 1);
for k = 1:numel(checkpointTimes)
    [~, idx(k)] = min(abs(tout - checkpointTimes(k)));
end
timeOut = tout(idx).';
deltaOut = xall(idx, 1:ngen);
omegaAbs = xall(idx, 16:(15+ngen));
omegaOut = omegaAbs - (omegaAbs * preset.m / sum(preset.m)) * ones(1, ngen);
thetaOut = xall(idx, 4:9);
voltageOut = xall(idx, 10:15);
residualOut = zeros(numel(idx), 1);
for k = 1:numel(idx)
    fval = F_3M9B_SP_DAE(xall(idx(k), :).', "fault1");
    residualOut(k) = norm([fval(4:8); fval(10:14)]);
end

ref = struct();
ref.name = 'spm_numerical_v2';
ref.schema_version = '1.0';
ref.status = 'AVAILABLE';
ref.reason = '';
ref.metadata = struct('case', 'case9_v2', 'fault', 'F9', ...
    'source_matlab_commit', '035f1475fd92e5639ff9b7fb78eb678ed2976e1c', ...
    'matlab_version', version('-release'), 'created_at', datestr(now, 30), ...
    'checkpoint_request', checkpointTimes, 'actual_time', timeOut, ...
    'tunit', Tunit, 'network_context', 'fault1', ...
    'network_placeholder_bus', faultbus, 'omega_frame', 'coi_relative');
ref.evidence = struct('checks_passed', numel(idx), 'checks_total', numel(idx), ...
    'max_error', max(residualOut), 'kind', 'matlab_native_fault1_dae', ...
    'note', '同一 fault1 DAE 的固定检查点；尚未声称与 Python 逐变量误差已通过');
ref.arrays = struct('time', timeOut, 'delta_gen', deltaOut, ...
    'omega_gen', omegaOut, 'theta_net', thetaOut, 'voltage_net', voltageOut, ...
    'algebraic_residual', residualOut, ...
    'phase', {repmat({'fault-on'}, 1, numel(idx))});

outfile = fullfile(projectRoot, '..', 'validation', 'references', 'spm_numerical_v2.json');
fid = fopen(outfile, 'w');
if fid < 0, error('SPM_FAULT_EXPORT_IO: cannot open %s', outfile); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, jsonencode(ref), 'char');
fprintf('EXPORT_SPM_FAULT_DAE_OK points=%d max_residual=%.6g file=%s\n', ...
    numel(idx), max(residualOut), outfile);
