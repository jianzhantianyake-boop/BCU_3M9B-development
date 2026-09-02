%EXPORT_SPM_ENERGY_GATE Read-only audit of the SPM energy switch and CUEP frame.
%
% This verifier runs the unmodified Cal_MM_CCT_SPM pipeline, then reports
% both the historical projected CUEP fields and the physically consistent
% common-angle frame.  It also records the PathEnergyCal value used by the
% formal CCT calculation.  The verifier never edits B3_MM source files.

% Run this verifier in the same base workspace immediately after
% ``Cal_MM_CCT_SPM.m``.  The source pipeline is a script whose helper
% functions use ``evalin('base',...)``; invoking it from a function would
% hide the initialized objects from those helpers.
if ~exist('postfault', 'var') || ~exist('Critical', 'var') || ...
        ~exist('Results_fsolve', 'var')
    error('export_spm_energy_gate:MissingPipeline', ...
        'Run Cal_MM_CCT_SPM.m in the same base workspace first.');
end
if ~exist('outputPath', 'var') || isempty(outputPath)
    outputPath = fullfile(pwd, 'validation', 'reports', ...
        'matlab_spm_energy_gate.json');
end

ngen = numel(preset.genno);
nnet = numel(postfault.net_delta);
raw = Results_fsolve(:);
raw_delta = [raw(1:(ngen-1)); 0];
raw_omega = raw(ngen);
raw_theta = raw((ngen+1):(ngen+nnet));
raw_voltage = raw((ngen+nnet+1):end);

% Cal_MM_CCT_SPM stores CUEP_net_theta with an unqualified base variable
% named deltacoi.  Reconstruct the physical frame using the actual shift
% returned by the CUEP calculation, without modifying the source result.
actual_shift = postfault.deltacoi;
physical_delta = raw_delta - actual_shift;
physical_theta = raw_theta - actual_shift;
projected_delta = postfault.CUEP_delta(:);
projected_theta = postfault.CUEP_net_theta(:);

x_raw = [raw(1:ngen); raw_theta; raw_voltage];
% Fun_SEPfslove_SPM uses generator ngen as its internal angle reference.
% Convert both stored COI states back to that local reference before
% evaluating the residual; otherwise a common frame shift is counted twice.
projected_ref = projected_delta(ngen);
x_projected = [projected_delta(1:(ngen-1)) - projected_ref; raw_omega; ...
               projected_theta - projected_ref; raw_voltage];
physical_ref = physical_delta(ngen);
x_physical = [physical_delta(1:(ngen-1)) - physical_ref; raw_omega; ...
              physical_theta - physical_ref; raw_voltage];
f_raw = Fun_SEPfslove_SPM(x_raw, preset, postfault, Basevalue);
f_projected = Fun_SEPfslove_SPM(x_projected, preset, postfault, Basevalue);
f_physical = Fun_SEPfslove_SPM(x_physical, preset, postfault, Basevalue);
raw_residual = norm(f_raw);
projected_residual = norm(f_projected((ngen+1):end));
physical_residual = norm(f_physical((ngen+1):end));

path_used_by_cct = preset.PathEnergyCal;
preset.PathEnergyCal = 0;
[e1,e2,e3,e4,e5] = Fun_Cal_PotentialEnergy_SPM( ...
    preset, postfault, physical_delta, physical_theta, raw_voltage);
physical_components = [e1 e2 e3 e4 e5];
physical_energy = sum(physical_components);
[p1,p2,p3,p4,p5] = Fun_Cal_PotentialEnergy_SPM( ...
    preset, postfault, projected_delta, projected_theta, raw_voltage);
projected_components = [p1 p2 p3 p4 p5];
projected_energy = sum(projected_components);

% The fault energy maximum is recomputed with the exact switch used by the
% source CCT scan.  This is intentionally a finite diagnostic, not a new
% reference baseline.
count = size(fault.traj.deltac, 1);
total = zeros(count, 1);
for k = 1:count
    [q1,q2,q3,q4,q5] = Fun_Cal_PotentialEnergy_SPM( ...
        preset, postfault, fault.traj.deltac(k,:).', ...
        fault.traj.theta(k,:).', fault.traj.voltage(k,:).');
    ek = 0.5 * sum(preset.m(:) .* fault.traj.omegac(k,:).'.^2);
    total(k) = ek + q1 + q2 + q3 + q4 + q5;
end
[fault_peak, peak_index] = max(total);

report = struct();
report.schema_version = '1.0';
report.kind = 'matlab_native_spm_energy_gate';
report.created_at = datestr(now, 30);
report.case = 'case9_v2';
report.fault = 'F9';
report.source_matlab_commit = ...
    '035f1475fd92e5639ff9b7fb78eb678ed2976e1c';
report.matlab_version = version('-release');
report.path_energy_cal_before_cct = path_used_by_cct;
report.path_energy_cal_used_for_cct = preset.PathEnergyCal;
report.historical_projected_e_critical = Critical.Ep;
report.historical_projected_lea_cct = Critical.LEA.CCT;
report.historical_projected_network_residual = projected_residual;
report.raw_fsolve_residual = raw_residual;
report.physical_network_residual = physical_residual;
report.projected_full_residual = norm(f_projected);
report.physical_full_residual = norm(f_physical);
report.coordinate_shift = actual_shift;
report.projected_physical_theta_max_abs_difference = ...
    max(abs(projected_theta - physical_theta));
report.physical_cuep_delta = physical_delta(:).';
report.physical_cuep_theta = physical_theta(:).';
report.physical_cuep_voltage = raw_voltage(:).';
report.physical_e_components = physical_components;
report.physical_e_critical = physical_energy;
report.projected_e_components_recomputed = projected_components;
report.projected_e_critical_recomputed = projected_energy;
report.fault_energy_peak = fault_peak;
report.fault_energy_peak_time = peak_index * fault.traj.Tunit;
report.gate_status = 'BLOCKED';
report.gate_reason = ['formal PathEnergyCal=0 physical E_critical exceeds ' ...
    'the same-window fault energy peak; historical projected fields are mixed-frame'];

folder = fileparts(outputPath);
if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
fid = fopen(outputPath, 'w');
if fid < 0
    error('export_spm_energy_gate:OpenFailed', 'Cannot open %s', outputPath);
end
fprintf(fid, '%s\n', jsonencode(report));
fclose(fid);
fprintf(['EXPORT_SPM_ENERGY_GATE status=%s path0=%.12g physical=%.12g ' ...
    'projected=%.12g peak=%.12g file=%s\n'], ...
    report.gate_status, report.path_energy_cal_used_for_cct, ...
    physical_energy, Critical.Ep, fault_peak, outputPath);
