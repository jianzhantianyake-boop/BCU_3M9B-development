%PROBE_SPM_CANDIDATE Read-only native check of a candidate SPM equilibrium.
%
% This script is deliberately separate from B3_MM.  It initializes the
% unmodified MATLAB SPM reference, solves Fun_SEPfslove_SPM from one explicit
% candidate, and exports the physically consistent COI-frame state and five
% potential-energy terms.  It does not overwrite any source result.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_bcu_paths(projectRoot);
Cal_MM_Static_SPM;

% Candidate supplied by the independent Python reduced type-1 search.  The
% network values are only an initial guess; MATLAB solves them again.
candidate_delta = [1.0325434376988165, -2.634499382702954, ...
                   -2.4942480047206974].';
candidate_theta = [1.0784924938902984, 1.488049206538747, ...
                   1.0130878072957497, -2.9304977459391135, ...
                   -2.866402046286545, -2.7283731405375797].';
candidate_voltage = [0.6115351306031335, 0.29128959201760224, ...
                     0.59290880685731, 0.5805594920170605, ...
                     0.6514770252735232, 0.7739313039940403].';

ngen = numel(preset.genno);
nbus = size(postfault.Yfull_mod, 1);
nnet = nbus - ngen;
candidate_ref = candidate_delta(ngen);
x_init = [candidate_delta(1:(ngen-1)) - candidate_ref; 0; ...
          candidate_theta - candidate_ref; candidate_voltage];
options = optimset('TolFun', 1e-12, 'MaxFunEvals', 1e5, ...
                   'MaxIter', 1e5, 'Display', 'off', 'TolX', 1e-10);
results = fsolve(@(x)Fun_SEPfslove_SPM(x, preset, postfault, Basevalue), ...
                 x_init, options);

raw_delta = [results(1:(ngen-1)); 0];
shift = raw_delta' * preset.m(:) / sum(preset.m(:));
physical_delta = raw_delta - shift;
physical_theta = results((ngen+1):(ngen+nnet)) - shift;
physical_voltage = results((ngen+nnet+1):end);
x_check = [physical_delta(1:(ngen-1)) - physical_delta(ngen); ...
           results(ngen); physical_theta - physical_delta(ngen); ...
           physical_voltage];
residual = norm(Fun_SEPfslove_SPM(x_check, preset, postfault, Basevalue));
preset.PathEnergyCal = 0;
[e1,e2,e3,e4,e5] = Fun_Cal_PotentialEnergy_SPM( ...
    preset, postfault, physical_delta, physical_theta, physical_voltage);

% Recompute the same-window fault energy curve used by the source CCT
% routine, rather than reusing a Python trajectory or a historical scalar.
Tfault = 0.5;
Tunit = 1e-4;
delta0 = prefault.SEP_delta;
omega0 = prefault.SEP_omegapu * Basevalue.omegab;
[fault_delta, fault_omega, fault_omegac, fault_theta, fault_voltage, ...
    ~, ~] = Fun_Cal_Exitpoint_SPM(Tfault, Tunit, fault, postfault, ...
    preset, delta0, omega0, prefault.net_delta, prefault.net_voltage, Basevalue);
fault_energy = nan(size(fault_delta, 1), 1);
for k = 1:size(fault_delta, 1)
    [q1,q2,q3,q4,q5] = Fun_Cal_PotentialEnergy_SPM( ...
        preset, postfault, fault_delta(k,:).', fault_theta(k,:).', ...
        fault_voltage(k,:).');
    omega_coi = fault_omegac(k,:).';
    fault_energy(k) = 0.5 * sum(preset.m(:) .* omega_coi.^2) + ...
        q1 + q2 + q3 + q4 + q5;
end
[fault_peak, peak_index] = max(fault_energy);
cross_index = find(fault_energy >= (e1 + e2 + e3 + e4 + e5), 1, 'first');
if isempty(cross_index)
    fault_cct = NaN;
    fault_cct_found = false;
else
    fault_cct = (cross_index - 1) * Tunit;
    fault_cct_found = true;
end

report = struct();
report.schema_version = '1.0';
report.kind = 'matlab_native_spm_candidate_probe';
report.created_at = datestr(now, 30);
report.case = 'case9_v2';
report.fault = 'F9';
report.source_matlab_commit = ...
    '035f1475fd92e5639ff9b7fb78eb678ed2976e1c';
report.matlab_version = version('-release');
report.raw_solver_residual = norm(Fun_SEPfslove_SPM(results, preset, postfault, Basevalue));
report.physical_solver_residual = residual;
report.coordinate_shift = shift;
report.cuep_delta = physical_delta(:).';
report.cuep_net_theta = physical_theta(:).';
report.cuep_net_voltage = physical_voltage(:).';
report.energy_components = [e1 e2 e3 e4 e5];
report.e_critical = sum(report.energy_components);
report.voltage_positive = all(physical_voltage > 1e-4);
report.fault_energy_points = numel(fault_energy);
report.fault_energy_peak = fault_peak;
report.fault_energy_peak_time = (peak_index - 1) * Tunit;
report.fault_energy_cct = fault_cct;
report.fault_energy_cct_found = fault_cct_found;

% Cal_MM_Static_SPM clears its own temporary workspace; recompute the script
% location instead of relying on the pre-initialization variable.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
outputPath = fullfile(projectRoot, '..', 'validation', 'reports', ...
    'matlab_spm_second_candidate.json');
folder = fileparts(outputPath);
if ~isfolder(folder), mkdir(folder); end
fid = fopen(outputPath, 'w');
if fid < 0, error('probe_spm_candidate:OpenFailed', 'Cannot open %s', outputPath); end
fprintf(fid, '%s\n', jsonencode(report));
fclose(fid);
fprintf(['PROBE_SPM_CANDIDATE residual=%.6g Ecritical=%.12g ' ...
         'positive=%d file=%s\n'], residual, report.e_critical, ...
        report.voltage_positive, outputPath);
