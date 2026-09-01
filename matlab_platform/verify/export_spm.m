% =========================================================================
% export_spm.m -- 导出 spm_cct 路径的 MATLAB 参考(SPM 能量法四件套结果 + 势能辅助量)
% 目标量: CUEP(发电机角+网络角+电压)、SEP 网络态、SPM 临界能量、CCT(LEA)、以及势能公式
% 所需辅助矩阵(Yfull_mod / Sload / Transform / genno) 供 Python 逐项对齐。
% =========================================================================
global VDIR
VDIR = fileparts(mfilename('fullpath'));
projectRoot = fileparts(VDIR);
addpath(projectRoot); setup_bcu_paths(projectRoot);
set(0, 'DefaultFigureVisible', 'off');
cd(projectRoot);

run(fullfile('B3_MM', 'Cal_MM_CCT_SPM.m'));   % 内部 clear + Cal_MM_Static_SPM + 全流程
global VDIR

ref = struct();
% --- CUEP (发电机角 + 网络角/电压) ---
ref.CUEP_delta       = postfault.CUEP_delta(:).';
ref.CUEP_net_theta   = postfault.CUEP_net_theta(:).';
ref.CUEP_net_voltage = postfault.CUEP_net_voltage(:).';
% --- SEP (发电机角 + 网络角/电压) ---
ref.SEP_delta        = postfault.SEP_delta(:).';
ref.SEP_net_delta    = postfault.net_delta(:).';
ref.SEP_net_voltage  = postfault.net_voltage(:).';
% --- 能量与 CCT ---
ref.CCT_LEA          = Critical.LEA.CCT;
ref.E_critical       = Critical.Ep;
ref.escape_deltac    = escape.deltac(:).';
ref.MGP_deltac       = MGP.detac_MGP(:).';
% --- 势能公式辅助量 ---
ref.Yfull_mod        = postfault.Yfull_mod;         % 重排后完整导纳(G/B)
ref.Transform        = postfault.Transform(:).';    % 母线重排映射
ref.genno            = preset.genno(:).';
ref.Sload            = preset.Sload;                % 负荷 [busno, P, Q]
ref.Pmpu             = preset.Pmpu(:).';
ref.Epu              = preset.Epu(:).';
ref.ngen             = preset.ngen;
ref.nbus             = preset.nbus;

% --- 故障轨迹(用于 CCT 逐点能量对比) ---
ref.traj_deltac  = fault.traj.deltac;
ref.traj_omegac  = fault.traj.omegac;
ref.traj_theta   = fault.traj.theta;     % 网络角(含故障母线补 0)
ref.traj_voltage = fault.traj.voltage;   % 网络电压
ref.traj_Tunit   = fault.traj.Tunit;

save(fullfile(VDIR, 'baseline_spm.mat'), 'ref', '-v7');
fprintf('EXPORT_SPM_OK CCT_LEA=%.4f E_crit=%.4f nbus=%d ngen=%d\n', ...
        ref.CCT_LEA, ref.E_critical, ref.nbus, ref.ngen);
