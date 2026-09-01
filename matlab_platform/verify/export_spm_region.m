% =========================================================================
% export_spm_region.m -- 导出 SPM 平衡点集合的紧凑 MATLAB 参考
%
% 只读调用原始 Statable_Region_SPM.m；不修改 B3_MM 中的方程、参数或判据。
% 输出 validation/references/spm_region_v1.json 所需的平衡点状态、分类和
% 分支标识。导出本身不是 Python/MATLAB 交叉验证，evidence.checks_* 保持 0。
% =========================================================================
thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisDir);
addpath(projectRoot); setup_bcu_paths(projectRoot);
set(0, 'DefaultFigureVisible', 'off');
cd(projectRoot);

run(fullfile('B3_MM', 'Cal_MM_Static_SPM.m'));
run(fullfile('B3_MM', 'Statable_Region_SPM.m'));

% ep_set 由原始脚本写入 base workspace。保留有限、正电压的记录。
if ~exist('ep_set', 'var') || isempty(ep_set)
    error('SPM_REGION_EXPORT_EMPTY: Statable_Region_SPM produced no equilibria');
end

m = preset.m(:);
ngen = preset.ngen;
nbus = preset.nbus;
nnet = nbus - ngen;
keep = false(1, length(ep_set));
for k = 1:length(ep_set)
    keep(k) = isfield(ep_set(k), 'xep') && isfield(ep_set(k), 'delta_net_ep') && ...
              isfield(ep_set(k), 'voltage_net_ep') && ...
              numel(ep_set(k).xep) == ngen-1 && ...
              numel(ep_set(k).delta_net_ep) == nnet && ...
              numel(ep_set(k).voltage_net_ep) == nnet && ...
              all(isfinite(ep_set(k).voltage_net_ep(:))) && ...
              all(ep_set(k).voltage_net_ep(:) > 1e-4);
end
ep_set = ep_set(keep);
nep = length(ep_set);
if nep == 0
    error('SPM_REGION_EXPORT_INVALID: no finite positive-voltage equilibria');
end

delta_gen = zeros(nep, ngen);
theta_net = zeros(nep, nnet);
voltage_net = zeros(nep, nnet);
flag = zeros(nep, 1);
flag_reduce = zeros(nep, 1);
equilibrium_type = cell(nep, 1);
branch_id = cell(nep, 1);
for k = 1:nep
    xep = ep_set(k).xep(:);
    d1 = -(m(2:ngen)' * xep) / m(1);
    d = [d1; xep];
    d = d - (m' * d) / sum(m);       % 明确记录为 COI 坐标
    delta_gen(k, :) = d(:).';
    theta_net(k, :) = ep_set(k).delta_net_ep(:).';
    voltage_net(k, :) = ep_set(k).voltage_net_ep(:).';
    if isfield(ep_set(k), 'flag'), flag(k) = ep_set(k).flag; end
    if isfield(ep_set(k), 'flag_reduce'), flag_reduce(k) = ep_set(k).flag_reduce; end
    if flag_reduce(k) == 0
        equilibrium_type{k} = 'SEP';
    elseif flag_reduce(k) == 1
        equilibrium_type{k} = 'type-1';
    else
        equilibrium_type{k} = sprintf('type-%d', flag_reduce(k));
    end
    branch_id{k} = sprintf('matlab-spm-%04d', k);
end

ref = struct();
ref.name = 'spm_region_v1';
ref.schema_version = '1.0';
ref.status = 'AVAILABLE';
ref.reason = 'MATLAB 原生 Statable_Region_SPM 平衡点导出；尚未与 Python 逐点对照';
ref.metadata = struct('case', 'case9_v2', 'fault', 'F9', ...
    'source_matlab_commit', '035f1475fd92e5639ff9b7fb78eb678ed2976e1c', ...
    'matlab_version', version('-release'), 'created_at', datestr(now, 30));
ref.evidence = struct('checks_passed', 0, 'checks_total', 0, ...
    'kind', 'matlab_native_export', ...
    'note', '导出状态不等于 MATLAB/Python 交叉验证');
ref.arrays = struct('delta_gen', delta_gen, 'theta_net', theta_net, ...
    'voltage_net', voltage_net, 'flag', flag, 'flag_reduce', flag_reduce, ...
    'equilibrium_type', {equilibrium_type}, 'branch_id', {branch_id});

% 原始 Statable_Region_SPM.m 会 clear 工作区变量，因此在返回后重新解析路径。
thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisDir);
outfile = fullfile(projectRoot, '..', 'validation', 'references', 'spm_region_v1.json');
jsonText = jsonencode(ref);
fid = fopen(outfile, 'w');
if fid < 0, error('SPM_REGION_EXPORT_IO: cannot open output'); end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, jsonText, 'char');
fprintf('EXPORT_SPM_REGION_OK nep=%d file=%s\n', nep, outfile);
