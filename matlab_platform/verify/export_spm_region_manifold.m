% =========================================================================
% export_spm_region_manifold.m -- 导出 MATLAB SPM 稳定流形固定检查点
%
% 只读复用 Statable_Region_SPM 的平衡点搜索和 f_reducedstate_SPM_backward，
% 将 type-1 UEP 的正/负稳定子空间扰动轨迹压缩为固定检查点。该文件不修改
% B3_MM 源码；导出本身不是跨平台通过证据，Python 对照由 run_full_xval.py 完成。
% =========================================================================
thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisDir);
addpath(projectRoot); setup_bcu_paths(projectRoot);
set(0, 'DefaultFigureVisible', 'off');
cd(projectRoot);

run(fullfile('B3_MM', 'Cal_MM_Static_SPM.m'));
run(fullfile('B3_MM', 'Statable_Region_SPM.m'));

% Statable_Region_SPM starts with ``clear`` and therefore removes the
% exporter-local path variables.  Re-resolve them before writing the output.
thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisDir);
addpath(projectRoot);
cd(projectRoot);

if ~exist('ep_set', 'var') || isempty(ep_set)
    error('SPM_REGION_MANIFOLD_EMPTY: no equilibria from Statable_Region_SPM');
end

m = preset.m(:);
ngen = preset.ngen;
nbus = preset.nbus;
nnet = nbus - ngen;
sampleTimes = [0, 0.25, 0.5, 0.75, 1.0];
perturb = 1e-2;

type1 = find(arrayfun(@(s) isfield(s, 'flag_reduce') && s.flag_reduce == 1, ep_set));
if isempty(type1)
    error('SPM_REGION_MANIFOLD_NO_TYPE1: no type-1 UEP');
end

branchCount = 2 * numel(type1);
deltaOut = zeros(branchCount, numel(sampleTimes), ngen);
thetaOut = zeros(branchCount, numel(sampleTimes), nnet);
voltageOut = zeros(branchCount, numel(sampleTimes), nnet);
branchSigns = zeros(branchCount, 1);
branchIds = cell(branchCount, 1);
residualOut = zeros(branchCount, numel(sampleTimes));
uepDelta = zeros(numel(type1), ngen);
uepTheta = zeros(numel(type1), nnet);
uepVoltage = zeros(numel(type1), nnet);
perturbVectors = zeros(numel(type1), 2);

row = 0;
options = odeset('Mass', diag([ones(2,1); ones(2*nnet,1)*1e-15]), ...
    'RelTol', 1e-10, 'AbsTol', [1e-8*ones(1,2), 1e-12*ones(1,2*nnet)]);
for q = 1:numel(type1)
    s = ep_set(type1(q));
    xep = s.xep(:);
    v = real(s.v_reduce(:,1));
    if norm(v) == 0 || any(~isfinite(v))
        error('SPM_REGION_MANIFOLD_BAD_VECTOR');
    end
    v = v / norm(v);
    d1 = -(m(2:ngen)' * xep) / m(1);
    d = [d1; xep];
    d = d - (m' * d) / sum(m);
    uepDelta(q,:) = d(:).';
    uepTheta(q,:) = s.delta_net_ep(:).';
    uepVoltage(q,:) = s.voltage_net_ep(:).';
    perturbVectors(q,:) = v(:).';

    for signValue = [-1, 1]
        row = row + 1;
        xp = xep + signValue * perturb * v;
        d1p = -(m(2:ngen)' * xp) / m(1);
        deltacc = [d1p; xp];
        net0 = [s.delta_net_ep(:); s.voltage_net_ep(:)];
        [netInit, ~, exitflag] = fsolve(@(z) Fun_AEfslove_SPM(z, deltacc, preset, "postfault"), ...
            net0, optimset('TolFun',1e-12,'TolX',1e-10,'Display','off'));
        if exitflag <= 0 || any(~isfinite(netInit))
            error('SPM_REGION_MANIFOLD_INIT: branch=%d fsolve failed', row);
        end
        [~, xall] = ode15s(@(t,x) f_reducedstate_SPM_backward(x), sampleTimes, ...
            [xp; netInit], options);
        branchSigns(row) = signValue;
        branchIds{row} = sprintf('matlab-spm-stable-%04d-%+d', q, signValue);
        for k = 1:numel(sampleTimes)
            x = xall(k,:).';
            d2 = x(1:2);
            d1x = -(m(2:ngen)' * d2) / m(1);
            dx = [d1x; d2];
            dx = dx - (m' * dx) / sum(m);
            deltaOut(row,k,:) = dx(:);
            thetaOut(row,k,:) = x(3:(2+nnet));
            voltageOut(row,k,:) = x((3+nnet):(2+2*nnet));
            % Only the algebraic P/Q constraint is a residual for this
            % reference.  The full backward vector also contains the
            % (deliberately nonzero) manifold tangent and must not be used as
            % a convergence metric.
            residualOut(row,k) = norm(Fun_AEfslove_SPM(x((3):(2+2*nnet)), ...
                dx, preset, "postfault"));
        end
    end
end

ref = struct();
ref.name = 'spm_region_manifold_v1';
ref.schema_version = '1.0';
ref.status = 'AVAILABLE';
ref.reason = 'MATLAB 原生稳定流形固定检查点；尚未与 Python 逐点对照';
ref.metadata = struct('case', 'case9_v2', 'fault', 'F9', ...
    'source_matlab_commit', '035f1475fd92e5639ff9b7fb78eb678ed2976e1c', ...
    'matlab_version', version('-release'), 'created_at', datestr(now, 30), ...
    'sample_times', sampleTimes, 'perturb', perturb, 'duration', sampleTimes(end));
ref.evidence = struct('checks_passed', 0, 'checks_total', 0, ...
    'kind', 'matlab_native_export', 'max_error', max(residualOut(:)), ...
    'note', '导出状态不等于 MATLAB/Python 交叉验证');
ref.arrays = struct('sample_time', sampleTimes, 'branch_sign', branchSigns, ...
    'branch_id', {branchIds}, 'uep_delta_gen', uepDelta, ...
    'uep_theta_net', uepTheta, 'uep_voltage_net', uepVoltage, ...
    'perturb_vectors', perturbVectors, 'delta_gen', deltaOut, ...
    'theta_net', thetaOut, 'voltage_net', voltageOut, ...
    'residual_norm', residualOut);

outfile = fullfile(projectRoot, '..', 'validation', 'references', 'spm_region_manifold_v1.json');
fid = fopen(outfile, 'w');
if fid < 0, error('SPM_REGION_MANIFOLD_IO: cannot open output'); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, jsonencode(ref), 'char');
fprintf('EXPORT_SPM_REGION_MANIFOLD_OK branches=%d points=%d max_residual=%.6g file=%s\n', ...
    branchCount, numel(sampleTimes), max(residualOut(:)), outfile);
