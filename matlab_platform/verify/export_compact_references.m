function manifest = export_compact_references(outputRoot)
%EXPORT_COMPACT_REFERENCES 导出不包含大型工作区的紧凑 JSON 参考。
% 该工具只写 validation/references 下的 JSON，不修改 B3_MM 方程，也不提交
% MATPOWER 或 .mat。当前 MATLAB 环境若无法运行某条管线，会生成 BLOCKED/
% UNVERIFIED 记录，而不是用零值填充。
if nargin < 1 || isempty(outputRoot)
    outputRoot = fullfile(fileparts(fileparts(mfilename('fullpath'))), '..', 'validation', 'references');
end
if ~isfolder(outputRoot), mkdir(outputRoot); end
sourceCommit = '035f1475fd92e5639ff9b7fb78eb678ed2976e1c';
manifest = struct('schema_version', '1.0', 'created_at', datestr(now, 30), ...
    'source_matlab_commit', sourceCommit, 'entries', []);
names = {'reduced_cct_v1','reduced_numerical_v1','reduced_region_v1', ...
    'two_machine_3d_v1','two_machine_gfl_v1','spm_cct_v1'};
for k = 1:numel(names)
    rec = struct('name', names{k}, 'schema_version', '1.0', ...
        'status', 'UNVERIFIED', 'reason', '尚未在当前 MATLAB 会话导出该路径的紧凑参考', ...
        'metadata', struct('case', 'case9_v2', 'fault', 'F9', ...
        'source_matlab_commit', sourceCommit, 'created_at', datestr(now, 30)), ...
        'arrays', struct());
    write_json(fullfile(outputRoot, [names{k} '.json']), rec);
    manifest.entries = [manifest.entries; struct('name', names{k}, ...
        'path', [names{k} '.json'], 'status', rec.status)]; %#ok<AGROW>
end
spm = struct('name', 'spm_cct_v1', 'schema_version', '1.0', ...
    'status', 'BLOCKED', 'reason', '未在本次会话运行 Cal_MM_CCT_SPM', ...
    'metadata', struct('case', 'case9_v2', 'fault', 'F9', ...
    'source_matlab_commit', sourceCommit, 'created_at', datestr(now, 30)), ...
    'arrays', struct());
try
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(projectRoot); setup_bcu_paths(projectRoot);
    if evalin('base', 'exist(''postfault'', ''var'') && exist(''Critical'', ''var'')')
        pf = evalin('base', 'postfault'); cr = evalin('base', 'Critical');
        spm.status = 'AVAILABLE'; spm.reason = '';
        spm.arrays = struct('sep_delta', pf.SEP_delta(:).', 'cuep_delta', pf.CUEP_delta(:).', ...
            'sep_net_theta', pf.net_delta(:).', 'sep_net_voltage', pf.net_voltage(:).', ...
            'cuep_net_theta', pf.CUEP_net_theta(:).', 'cuep_net_voltage', pf.CUEP_net_voltage(:).', ...
            'e_critical', cr.Ep, 'lea_cct', cr.LEA.CCT);
    end
catch err
    spm.status = 'BLOCKED'; spm.reason = ['MATLAB SPM 导出失败: ' err.message];
end
write_json(fullfile(outputRoot, 'spm_cct_v1.json'), spm);
manifest.entries(end) = struct('name', 'spm_cct_v1', 'path', 'spm_cct_v1.json', 'status', spm.status);
write_json(fullfile(outputRoot, 'reference_manifest.json'), manifest);
disp(['COMPACT_REFERENCE_MANIFEST status=' spm.status ' output=' outputRoot]);
end

function write_json(path, value)
fid = fopen(path, 'w');
if fid < 0, error('BCU:ReferenceWrite', '无法写入 %s', path); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, jsonencode(value));
end
