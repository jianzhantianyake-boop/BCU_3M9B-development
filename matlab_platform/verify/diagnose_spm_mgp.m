function report = diagnose_spm_mgp(sourceFile, outputFile)
%DIAGNOSE_SPM_MGP Read-only diagnostics for the MATLAB SPM MGP/CUEP stage.
%
% This function is deliberately separate from B3_MM.  It observes the base
% workspace after Cal_MM_CCT_SPM has run and records the suspicious historical
% assignment in Fun_Cal_MGP_SPM.m without changing that implementation.  A
% missing base-workspace variable is recorded as unavailable rather than being
% replaced by a fabricated value.

if nargin < 1 || isempty(sourceFile)
    sourceFile = fullfile(fileparts(mfilename('fullpath')), '..', 'B3_MM', ...
        'Fun_Cal_MGP_SPM.m');
end
if nargin < 2 || isempty(outputFile)
    outputFile = fullfile(pwd, 'spm_mgp_diagnostic.json');
end

report = struct();
report.schema_version = '1.0';
report.created_at = datestr(now, 30);
report.source_file = local_portable_path(sourceFile);
report.source_file_exists = isfile(sourceFile);
report.suspected_copy_line_present = false;
if report.source_file_exists
    sourceText = fileread(sourceFile);
    report.suspected_copy_line_present = ~isempty(strfind( ...
        sourceText, 'voltage_MGP= theta_iter(no_MGP,:)'));
end

report.branch_executed = evalin('base', 'exist(''MGP'', ''var'') == 1');
report.mgp = local_field_snapshot('MGP', {'detac_MGP','theta_MGP', ...
    'voltage_MGP','num_Traj','flag_MGP'});
report.cuep = local_field_snapshot('postfault', {'CUEP_delta', ...
    'CUEP_omegapu','CUEP_net_theta','CUEP_net_voltage','CUEP_Perr'});
report.fsolve = local_workspace_snapshot({'x_init','Results_fsolve'});
report.physical_cuep = local_physical_cuep_snapshot();
report.notes = { ...
    '只读诊断：不修改 B3_MM 或 Fun_Cal_MGP_SPM.m。', ...
    '疑似复制错误行只记录是否存在，不据此直接认定 MATLAB 缺陷。', ...
    '若变量未保留在 base workspace，则对应字段为 unavailable。', ...
    'physical_cuep 使用 Results_fsolve 的原始网络角减去实际 CUEP COI 平移，单独评估正确 SPM 网络残差与势能。'};

parent = fileparts(outputFile);
if ~isempty(parent) && ~isfolder(parent)
    mkdir(parent);
end

function snapshot = local_physical_cuep_snapshot()
% Reconstruct the CUEP network state in one consistent COI frame and check
% it against the physical SPM P/Q equations.  This is diagnostic-only: it
% never changes the MATLAB B3_MM solver or its base-workspace variables.
snapshot = struct('available', false);
try
    hasVars = evalin('base', ['exist(''postfault'', ''var'') == 1 && ' ...
        'exist(''preset'', ''var'') == 1 && exist(''Results_fsolve'', ''var'') == 1']);
    if ~hasVars
        snapshot.unavailable_reason = 'postfault/preset/Results_fsolve unavailable';
        return;
    end
    postfault = evalin('base', 'postfault');
    preset = evalin('base', 'preset');
    raw = evalin('base', 'Results_fsolve');
    ngen = size(preset.genno, 1);
    nbus = size(postfault.Yfull_mod, 1);
    if ~isfield(postfault, 'CUEP_delta') || ~isfield(postfault, 'deltacoi') || ...
            ~isfield(postfault, 'CUEP_net_voltage') || numel(raw) < 2 * nbus - ngen
        snapshot.unavailable_reason = 'CUEP fields or raw fsolve width unavailable';
        return;
    end
    delta = postfault.CUEP_delta(:);
    raw_theta = raw((ngen + 1):nbus);
    voltage = postfault.CUEP_net_voltage(:);
    frame_shift = postfault.deltacoi;
    corrected_theta = raw_theta - frame_shift;
    Y = postfault.Yfull_mod;
    all_theta = [delta; corrected_theta];
    all_voltage = [preset.Epu(:); voltage];
    phasor = all_voltage .* exp(1i * all_theta);
    injection = phasor .* conj(Y * phasor);
    network_residual = [real(injection((ngen + 1):nbus)); ...
        imag(injection((ngen + 1):nbus))];
    [Ep1, Ep2, Ep3, Ep4, Ep5] = Fun_Cal_PotentialEnergy_SPM(...
        preset, postfault, delta, corrected_theta, voltage);
    Ep = [Ep1, Ep2, Ep3, Ep4, Ep5];
    snapshot.available = true;
    snapshot.frame_shift = frame_shift;
    snapshot.raw_net_theta = raw_theta;
    snapshot.corrected_net_theta = corrected_theta;
    snapshot.cuep_delta = delta;
    snapshot.cuep_net_voltage = voltage;
    snapshot.network_residual = norm(network_residual);
    snapshot.energy_components = Ep(:).';
    snapshot.e_critical = sum(Ep);
catch err
    snapshot.available = false;
    snapshot.unavailable_reason = err.message;
end
end
fid = fopen(outputFile, 'w');
if fid < 0
    error('diagnose_spm_mgp:OpenFailed', '无法写入诊断文件: %s', outputFile);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', jsonencode(report));
end

function snapshot = local_field_snapshot(varName, fields)
snapshot = struct();
existsVar = evalin('base', sprintf('exist(''%s'', ''var'') == 1', varName));
snapshot.available = existsVar;
if ~existsVar
    snapshot.unavailable_reason = 'base workspace 中不存在变量';
    return;
end
for k = 1:numel(fields)
    fieldName = fields{k};
    expression = sprintf('isstruct(%s) && isfield(%s, ''%s'')', ...
        varName, varName, fieldName);
    hasField = evalin('base', expression);
    if hasField
        value = evalin('base', sprintf('%s.%s', varName, fieldName));
        snapshot.(fieldName) = local_json_value(value);
    else
        snapshot.([fieldName '_available']) = false;
    end
end
end

function snapshot = local_workspace_snapshot(varNames)
snapshot = struct();
for k = 1:numel(varNames)
    name = varNames{k};
    existsVar = evalin('base', sprintf('exist(''%s'', ''var'') == 1', name));
    snapshot.([name '_available']) = existsVar;
    if existsVar
        snapshot.(name) = local_json_value(evalin('base', name));
    end
end
end

function value = local_json_value(raw)
% Keep diagnostics compact while preserving numeric vectors/matrices.
if isnumeric(raw) || islogical(raw)
    value = raw;
elseif ischar(raw) || isstring(raw)
    value = char(raw);
else
    value = sprintf('<%s unavailable for compact JSON>', class(raw));
end
end

function portable = local_portable_path(pathName)
% Avoid embedding a user's absolute Windows profile path in committed data.
portable = char(pathName);
marker = [filesep 'matlab_platform' filesep];
idx = strfind(portable, marker);
if ~isempty(idx)
    portable = portable(idx(end) + 1:end);
end
portable = strrep(portable, '\', '/');
end
