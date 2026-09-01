function record = export_spm_trajectory_checkpoints(outputFile, checkpointTimes)
%EXPORT_SPM_TRAJECTORY_CHECKPOINTS Export compact MATLAB SPM fault checkpoints.
%
% The function observes fault.traj created by Cal_MM_CCT_SPM and writes only
% fixed checkpoints, not the full long trajectory.  It is intentionally a
% verify/ utility and does not alter the B3_MM simulation.
if nargin < 1 || isempty(outputFile)
    outputFile = fullfile(pwd, 'spm_numerical_v1.json');
end
if nargin < 2 || isempty(checkpointTimes)
    checkpointTimes = [0, 0.1, 0.2, 0.2053, 0.21, 0.25, 0.3];
end

record = struct('name', 'spm_numerical_v1', 'schema_version', '1.0', ...
    'status', 'BLOCKED', 'reason', 'base workspace 中没有 fault.traj', ...
    'metadata', struct('case', 'case9_v2', 'fault', 'F9', ...
    'source_matlab_commit', '035f1475fd92e5639ff9b7fb78eb678ed2976e1c', ...
    'created_at', datestr(now, 30)), 'arrays', struct());
try
    if ~evalin('base', 'exist(''fault'', ''var'') && isfield(fault, ''traj'')')
        write_json(outputFile, record); return;
    end
    fault = evalin('base', 'fault');
    preset = evalin('base', 'preset');
    prefault = evalin('base', 'prefault');
    traj = fault.traj;
    n = size(traj.deltac, 1);
    t = (1:n).' * traj.Tunit;
    % Include a true t=0 record from the prefault SEP; subsequent entries are
    % nearest available points in fault.traj (never interpolated or zero-filled).
    idx = zeros(numel(checkpointTimes), 1);
    for k = 1:numel(checkpointTimes)
        if checkpointTimes(k) <= 0
            idx(k) = 0;
        else
            [~, idx(k)] = min(abs(t - checkpointTimes(k)));
        end
    end
    timeOut = zeros(numel(idx), 1);
    deltaOut = cell(numel(idx), 1);
    omegaOut = cell(numel(idx), 1);
    thetaOut = cell(numel(idx), 1);
    voltageOut = cell(numel(idx), 1);
    for k = 1:numel(idx)
        if idx(k) == 0
            timeOut(k) = 0;
            deltaOut{k} = prefault.SEP_delta(:).';
            omegaOut{k} = zeros(1, numel(prefault.SEP_delta));
            thetaOut{k} = prefault.net_delta(:).';
            voltageOut{k} = prefault.net_voltage(:).';
        else
            timeOut(k) = t(idx(k));
            deltaOut{k} = traj.deltac(idx(k), :);
            omegaOut{k} = traj.omegac(idx(k), :);
            thetaOut{k} = traj.theta(idx(k), :);
            voltageOut{k} = traj.voltage(idx(k), :);
        end
    end
    record.status = 'AVAILABLE';
    record.reason = '';
    record.metadata.checkpoint_request = checkpointTimes(:).';
    record.metadata.actual_time = timeOut(:).';
    record.metadata.tunit = traj.Tunit;
    record.metadata.escape_index = evalin('base', 'escape.tm');
    record.arrays = struct('time', timeOut(:).', 'delta_gen', {deltaOut}, ...
        'omega_gen', {omegaOut}, 'theta_net', {thetaOut}, ...
        'voltage_net', {voltageOut}, 'phase', 'fault-on');
catch err
    record.status = 'BLOCKED';
    record.reason = ['SPM trajectory export failed: ' err.message];
end
write_json(outputFile, record);
end

function write_json(path, value)
parent = fileparts(path);
if ~isempty(parent) && ~isfolder(parent), mkdir(parent); end
fid = fopen(path, 'w');
if fid < 0, error('BCU:ReferenceWrite', '无法写入 %s', path); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, jsonencode(value));
end
