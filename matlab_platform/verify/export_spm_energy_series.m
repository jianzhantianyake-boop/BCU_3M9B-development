%EXPORT_SPM_ENERGY_SERIES Read-only MATLAB SPM fault energy diagnostic.
%
% This verifier follows the native Cal_MM_CCT_SPM sequence up to the
% fault-on/postfault algebraic trajectory, then evaluates the exact
% Fun_Cal_CCT_Energy_SPM functional at every point.  It does not edit B3_MM
% and is intentionally separate from the compact reference exporters.

function report = export_spm_energy_series(outputPath)
if nargin < 1 || isempty(outputPath)
    outputPath = fullfile(pwd, 'validation', 'reports', 'matlab_spm_energy_series.json');
end
% Cal_MM_Static_SPM starts with ``clear`` and can clear variables in the
% workspace in which this verifier is run.  Keep the requested output path
% in base workspace so it survives that source script.
assignin('base', 'spm_energy_output_path', outputPath);
set(0, 'DefaultFigureVisible', 'off');
thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisDir);
addpath(projectRoot);
setup_bcu_paths(projectRoot);
cd(projectRoot);

% Cal_MM_Static_SPM is a source script and clears its caller workspace.
run(fullfile('B3_MM', 'Cal_MM_Static_SPM.m'));
outputPath = evalin('base', 'spm_energy_output_path');
% The legacy SPM functions read these objects via evalin('base',...).
% Publish only the initialized objects required by this read-only verifier.
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
[delta, omega, omegac, theta, voltage, exit_tm, dotproduct] = ...
    Fun_Cal_Exitpoint_SPM(Tfault, Tunit, fault, postfault, preset, ...
                          delta0, omega0, delta_net0, voltage_net0, Basevalue);

% Cal_MM_CCT_SPM resets this setting before its critical-energy scan.  Do the
% same here so this diagnostic is a strict comparison of the CCT functional,
% not the 20-segment MGP ray integral used during static initialization.
preset.PathEnergyCal = 0;

count = size(delta, 1);
kinetic = zeros(count, 1);
potential = zeros(count, 1);
total = zeros(count, 1);
components = zeros(count, 5);
for k = 1:count
    [e1, e2, e3, e4, e5] = Fun_Cal_PotentialEnergy_SPM( ...
        preset, postfault, delta(k, :).', theta(k, :).', voltage(k, :).');
    components(k, :) = [e1, e2, e3, e4, e5];
    kinetic(k) = 0.5 * sum(preset.m(:) .* omegac(k, :).'.^2);
    potential(k) = sum(components(k, :));
    total(k) = kinetic(k) + potential(k);
end

[peak, peakIndex] = max(total);
report = struct();
report.schema_version = '1.0';
report.kind = 'matlab_native_spm_fault_energy_series';
report.created_at = datestr(now, 30);
report.case = 'case9_v2';
report.fault = 'F9';
report.source_matlab_commit = '035f1475fd92e5639ff9b7fb78eb678ed2976e1c';
report.matlab_version = version('-release');
report.tfault = Tfault;
report.tunit = Tunit;
report.count = count;
report.exit_tm = exit_tm;
report.exit_time = exit_tm * Tunit;
report.peak_index = peakIndex;
report.peak_time = peakIndex * Tunit;
report.peak_energy = peak;
report.peak_kinetic = kinetic(peakIndex);
report.peak_potential = potential(peakIndex);
report.peak_components = components(peakIndex, :);
report.dotproduct_at_peak = dotproduct(peakIndex);
report.first_delta = delta(1, :);
report.first_theta = theta(1, :);
report.first_voltage = voltage(1, :);
report.first_components = components(1, :);
report.last_delta = delta(end, :);
report.last_theta = theta(end, :);
report.last_voltage = voltage(end, :);
report.last_components = components(end, :);
report.peak_delta = delta(peakIndex, :);
report.peak_theta = theta(peakIndex, :);
report.peak_voltage = voltage(peakIndex, :);
report.arrays = struct('time', (1:count) * Tunit, ...
    'kinetic', kinetic.', 'potential', potential.', 'total', total.');

folder = fileparts(outputPath);
if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
fid = fopen(outputPath, 'w');
if fid < 0
    error('export_spm_energy_series:OpenFailed', 'Cannot open %s', outputPath);
end
fprintf(fid, '%s\n', jsonencode(report));
fclose(fid);
fprintf('EXPORT_SPM_ENERGY_SERIES_OK count=%d peak=%.12g t=%.7g file=%s\n', ...
    count, peak, report.peak_time, outputPath);
end
