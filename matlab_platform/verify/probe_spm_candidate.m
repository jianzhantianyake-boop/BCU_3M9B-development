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
