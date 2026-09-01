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
report.notes = { ...
    '只读诊断：不修改 B3_MM 或 Fun_Cal_MGP_SPM.m。', ...
    '疑似复制错误行只记录是否存在，不据此直接认定 MATLAB 缺陷。', ...
    '若变量未保留在 base workspace，则对应字段为 unavailable。'};

parent = fileparts(outputFile);
if ~isempty(parent) && ~isfolder(parent)
    mkdir(parent);
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
    portable = portable(idx(end) + numel(filesep):end);
    portable = ['matlab_platform' filesep portable];
end
portable = strrep(portable, '\', '/');
end
