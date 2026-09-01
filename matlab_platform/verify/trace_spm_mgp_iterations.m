%TRACE_SPM_MGP_ITERATIONS Read-only trace of the MATLAB SPM MGP outer loop.
% This verifier mirrors Cal_MM_CCT_SPM's call sequence without editing B3_MM.
% It records each trajectory's start/end point, norm sequence, and ray update
% so a Python implementation can be compared to actual branch execution.

function report = trace_spm_mgp_iterations(outputPath)
if nargin < 1 || isempty(outputPath)
    outputPath = fullfile(pwd, 'validation', 'reports', 'spm_mgp_iterations.json');
end
set(0, 'DefaultFigureVisible', 'off');
assignin('base', 'spm_mgp_trace_output_path', outputPath);
run(fullfile('matlab_platform', 'B3_MM', 'Cal_MM_Static_SPM.m'));
outputPath = evalin('base', 'spm_mgp_trace_output_path');
% The legacy SPM helpers use evalin('base',...).  Keep this verifier's
% workspace isolated while explicitly publishing only the initialized state
% objects they require; no source platform file is modified.
assignin('base', 'Basevalue', Basevalue);
assignin('base', 'preset', preset);
assignin('base', 'prefault', prefault);
assignin('base', 'fault', fault);
assignin('base', 'postfault', postfault);

Tfault = 0.5;
Tunit = 1e-4;
delta0 = prefault.SEP_delta;
omega0 = prefault.SEP_omegapu * Basevalue.omegab;
delta_net0 = prefault.net_delta;
voltage_net0 = prefault.net_voltage;
[deltac,omega,omegac,theta,voltage,exit_tm,Dotproduct] = ...
    Fun_Cal_Exitpoint_SPM(Tfault,Tunit,fault,postfault,preset,delta0,omega0,...
                          delta_net0,voltage_net0,Basevalue);
escape.deltac = deltac(exit_tm,:);
escape.theta = theta(exit_tm,:);
escape.voltage = voltage(exit_tm,:);

trace = struct('trajectory', {}, 'start_delta', {}, 'end_delta', {}, ...
               'norms', {}, 'norm_min', {}, 'flag_mgp', {}, ...
               'update_flag', {}, 'update_delta', {}, 'update_energy', {}, ...
               'ray_max_energy', {}, 'ray_max_distance', {}, ...
               'ray_last_energy', {}, 'ray_sample_energy', {}, ...
               'ray_sample_components', {}, ...
               'ray_sample_indices', {}, 'ray_sample_vmin', {}, ...
               'ray_sample_theta', {}, 'ray_sample_voltage', {});
deltac_start = escape.deltac';
theta_start = escape.theta';
voltage_start = escape.voltage';
flag_MGP = 0;
traj_no = 0;
norm_min = 0;
while flag_MGP == 0 && traj_no < 1000
    [theta_start,voltage_start,flag_iter,n_iter,err] = ...
        Fun_AEiteration_SPM(theta_start,voltage_start,deltac_start,preset,...
                            Basevalue,"postfault",1e4,1e-10);
    [deltac_iter,theta_iter,voltage_iter,Normp,no_MGP,flag_MGP] = ...
        Fun_Cal_MGP_singletraj_SPM(deltac_start,theta_start,voltage_start,...
                                   1e-3,10,1e-5,postfault.Yfull_mod,preset);
    traj_no = traj_no + 1;
    item = struct();
    item.trajectory = traj_no;
    item.start_delta = deltac_start(:)';
    item.end_delta = deltac_iter(end,:)';
    item.norms = Normp(:)';
    item.norm_min = min(Normp);
    item.flag_mgp = flag_MGP;
    item.update_flag = 0;
    item.update_delta = [];
    item.update_energy = NaN;
    if flag_MGP == 0
        deltac_last = deltac_iter(end,:)';
        theta_lastpoint = theta_iter(end,:)';
        voltage_lastpoint = voltage_iter(end,:)';
        if traj_no == 1
            % Independently sample the first ray before calling the legacy
            % updater, so a Python implementation can distinguish a genuine
            % local maximum from a branch/path discrepancy.
            sep_delta = postfault.SEP_delta(:);
            ray = deltac_last - sep_delta;
            ray = ray / norm(ray);
            n_ray = fix(2*norm(deltac_last-sep_delta)/1e-3);
            ray_energy = NaN(n_ray+1,1);
            ray_components = NaN(n_ray+1,5);
            sample_idx = [1 501 1001 1801 2501 3501 5001];
            sample_idx = sample_idx(sample_idx <= n_ray+1);
            ray_point = sep_delta;
            ray_theta = postfault.net_delta(:);
            ray_voltage = postfault.net_voltage(:);
            ray_theta_samples = NaN(numel(sample_idx), numel(ray_theta));
            ray_voltage_samples = NaN(numel(sample_idx), numel(ray_voltage));
            [e1,e2,e3,e4,e5] = Fun_Cal_PotentialEnergy_SPM(...
                preset,postfault,ray_point,ray_theta,ray_voltage);
            ray_components(1,:) = [e1 e2 e3 e4 e5];
            ray_energy(1) = e1+e2+e3+e4+e5;
            sample_pos = find(sample_idx == 1, 1);
            if ~isempty(sample_pos)
                ray_theta_samples(sample_pos,:) = ray_theta(:)';
                ray_voltage_samples(sample_pos,:) = ray_voltage(:)';
            end
            for k = 1:n_ray
                ray_point = sep_delta + k*1e-3*ray;
                [ray_theta,ray_voltage,flag_ray,~,~] = Fun_AEiteration_SPM(...
                    ray_theta,ray_voltage,ray_point,preset,Basevalue,...
                    "postfault",1e4,1e-12);
                if flag_ray ~= 1
                    ray_energy = ray_energy(1:k);
                    break;
                end
                [e1,e2,e3,e4,e5] = Fun_Cal_PotentialEnergy_SPM(...
                    preset,postfault,ray_point,ray_theta,ray_voltage);
                ray_components(k+1,:) = [e1 e2 e3 e4 e5];
                ray_energy(k+1) = e1+e2+e3+e4+e5;
                sample_pos = find(sample_idx == k+1, 1);
                if ~isempty(sample_pos)
                    ray_theta_samples(sample_pos,:) = ray_theta(:)';
                    ray_voltage_samples(sample_pos,:) = ray_voltage(:)';
                end
            end
            [ray_max,ray_max_i] = max(ray_energy);
            item.ray_max_energy = ray_max;
            item.ray_max_distance = (ray_max_i-1)*1e-3;
            item.ray_last_energy = ray_energy(end);
            item.ray_sample_indices = sample_idx - 1;
            item.ray_sample_energy = ray_energy(sample_idx)';
            item.ray_sample_components = ray_components(sample_idx,:);
            item.ray_sample_vmin = NaN(size(sample_idx));
            item.ray_sample_theta = ray_theta_samples;
            item.ray_sample_voltage = ray_voltage_samples;
        else
            item.ray_max_energy = NaN;
            item.ray_max_distance = NaN;
            item.ray_last_energy = NaN;
            item.ray_sample_energy = [];
            item.ray_sample_components = [];
            item.ray_sample_indices = [];
            item.ray_sample_vmin = [];
            item.ray_sample_theta = [];
            item.ray_sample_voltage = [];
        end
        [deltac_update,theta_update,voltage_update,flag_update] = ...
            Fun_Cal_UpdateStartPoint_SPM(deltac_last,theta_lastpoint,...
                                         voltage_lastpoint,preset,postfault);
        item.update_flag = flag_update;
        item.update_delta = deltac_update(:)';
        if flag_update == 1
            [ep1,ep2,ep3,ep4,ep5] = Fun_Cal_PotentialEnergy_SPM(...
                preset,postfault,deltac_update,theta_update,voltage_update);
            item.update_energy = ep1+ep2+ep3+ep4+ep5;
        end
        deltac_starthis = deltac_start;
        deltac_start = deltac_update;
        theta_start = theta_update;
        voltage_start = voltage_update;
        if norm(deltac_start-deltac_starthis) < 1e-3
            flag_MGP = 1;
        elseif norm(deltac_start-deltac_starthis) > ...
                0.5*norm(deltac_starthis-postfault.SEP_delta)
            deltac_start = deltac_last;
            theta_start = theta_lastpoint;
            voltage_start = voltage_lastpoint;
        end
    end
    trace(end+1) = item; %#ok<AGROW>
end

report = struct();
report.exit_tm = exit_tm;
report.escape_delta = escape.deltac(:)';
report.trajectory_count = traj_no;
report.mgp_found = logical(flag_MGP);
report.mgp_delta = [];
report.mgp_theta = [];
report.mgp_voltage = [];
if ~isempty(trace)
    report.mgp_delta = trace(end).end_delta;
end
report.trace = trace;
folder = fileparts(outputPath);
if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
fid = fopen(outputPath, 'w');
if fid < 0, error('trace_spm_mgp_iterations:OpenFailed', 'Cannot open %s', outputPath); end
fprintf(fid, '%s\n', jsonencode(report));
fclose(fid);
end
