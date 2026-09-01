% =========================================================================
% run_bcu.m —— BCU_3M9B 仿真平台【主操作入口】（配置驱动版）
%
% 工作流：
%   1) 定位项目根、初始化路径
%   2) 读取 bcu_config()（你唯一需要编辑的参数文件）
%   3) 校验配置（维度/支路/步长…），有错即停并定位
%   4) 打印配置摘要，存参数快照到 results/（可追溯）
%   5) 按 cfg.mode 运行对应实验链路
%   6) 对 reduced_cct/numerical 自动做残差 + 能量法保守性自检
%
% 与旧入口的关系：run_bcu_beginner.m 仍可用（只改一行 mode）；本入口面向
% 需要微调物理/数值参数的场景，参数集中在 bcu_config.m，绝不改动核心方程。
% =========================================================================
clc;
thisFile = mfilename('fullpath');
projectRoot = fileparts(thisFile);
cd(projectRoot);
setup_bcu_paths(projectRoot);

bcu_override('clear');   % 单次入口不使用运行时覆盖（清除扫参可能残留的覆盖）
cfg = bcu_config();

% ----------------------------- 校验 -------------------------------------
fprintf('==================== BCU 配置校验 ====================\n');
[ok, msgs] = bcu_validate_config(cfg);
for k = 1:numel(msgs); fprintf('  %s\n', msgs{k}); end
if ~ok
    error('BCU:BadConfig', '配置校验未通过，请修正 bcu_config.m 后重试（见上方 [ERR] 行）。');
end

% --------------------------- 配置摘要 -----------------------------------
fprintf('\n==================== 本次实验配置 ====================\n');
fprintf('  实验模式 mode        : %s\n', cfg.mode);
fprintf('  案例/基频            : %s @ %g Hz\n', cfg.CaseName, cfg.f_base);
fprintf('  发电机数 ngen        : %d\n', numel(cfg.m));
fprintf('  惯性 m               : [%s]\n', num2str(cfg.m(:).', '%.4g '));
fprintf('  阻尼比 d/m           : [%s]\n', num2str(cfg.damping_ratio(:).', '%.4g '));
fprintf('  机械功率 Pm          : [%s]\n', num2str(cfg.Pm(:).', '%.4g '));
fprintf('  故障支路/位置        : [%d-%d] / pos %d\n', cfg.faultline(1), cfg.faultline(2), cfg.faultposition);
fprintf('  求解器/路径能量      : EquCal=%d, PathEnergyCal=%d\n', cfg.EquCal, cfg.PathEnergyCal);
fprintf('  积分 Tfault/Tunit    : %.3g s / %.1e s\n', cfg.Tfault, cfg.Tunit);

% --------------------------- 存参数快照 ---------------------------------
if cfg.save_snapshot
    resdir = fullfile(projectRoot, 'results');
    if ~isfolder(resdir); mkdir(resdir); end
    stamp = datestr(now, 'yyyymmdd_HHMMSS');
    save(fullfile(resdir, ['snapshot_' stamp '.mat']), 'cfg');
    fid = fopen(fullfile(resdir, ['snapshot_' stamp '.txt']), 'w');
    if fid > 0
        fn = fieldnames(cfg);
        fprintf(fid, '%% BCU 配置快照 %s\n', stamp);
        for k = 1:numel(fn)
            val = cfg.(fn{k});
            if isnumeric(val); vs = mat2str(val(:).'); else; vs = char(string(val)); end
            fprintf(fid, 'cfg.%s = %s;\n', fn{k}, vs);
        end
        fclose(fid);
    end
    fprintf('  参数快照已存        : results/snapshot_%s.{mat,txt}\n', stamp);
end
fprintf('=====================================================\n\n');

% --------------------------- 运行实验链路 -------------------------------
% 注意：Cal_MM_CCT*.m 会 clear 工作区并在内部重新调用 bcu_config()，
% 因此参数通过函数配置注入，不依赖此处的工作区变量。
switch cfg.mode
    case "reduced_cct"
        run(fullfile('B3_MM','Cal_MM_CCT.m'));
    case "reduced_numerical"
        run(fullfile('B3_MM','Cal_MM_CCT.m'));
        run(fullfile('B3_MM','NumSim_MM_Gridframe.m'));
    case "reduced_region"
        run(fullfile('B3_MM','Cal_MM_Static.m'));
        run(fullfile('B3_MM','Statable_Region.m'));
    case "spm_cct"
        run(fullfile('B3_MM','Cal_MM_CCT_SPM.m'));
    case "spm_numerical"
        run(fullfile('B3_MM','Cal_MM_CCT_SPM.m'));
        run(fullfile('B3_MM','NumSim_MM_Gridframe_SPM.m'));
    case "spm_region"
        run(fullfile('B3_MM','Cal_MM_Static_SPM.m'));
        run(fullfile('B3_MM','Statable_Region_SPM.m'));
    case "two_machine_region_3d"
        run(fullfile('B3_MM','Statable_Region_3D.m'));
    case "two_machine_region_3d_gfl"
        run(fullfile('B3_MM','Statable_Region_3D_GFL.m'));
end

% --------------------------- 运行后自检 ---------------------------------
% 工作区已被链路脚本 clear，故重新读取配置判断是否自检。
cfgEnd = bcu_config();
if cfgEnd.run_selfcheck && any(cfgEnd.mode == ["reduced_cct","reduced_numerical"])
    fprintf('\n============ 运行后自检（reduced CCT）============\n');
    try
        se = norm(postfault.SEP_Perr); ce = norm(postfault.CUEP_Perr);
        fprintf('  postfault SEP 残差   : %.2e  (%s)\n', se, passfail(se<1e-6));
        fprintf('  postfault CUEP 残差  : %.2e  (%s)\n', ce, passfail(ce<1e-6));
        lea = Critical.LEA.CCT; rea = Critical.REA.CCT;
        fprintf('  LEA-CCT / REA-CCT    : %.4f / %.4f s\n', lea, rea);
        fprintf('  能量法保守性 LEA<=REA: %s\n', passfail(lea <= rea + 1e-9));
    catch ME
        fprintf('  [自检跳过] %s\n', ME.message);
    end
    fprintf('=================================================\n');
end

fprintf('\nrun_bcu 完成。请检查命令窗口警告与图窗后再保存结果。\n');

function s = passfail(b)
    if b; s = 'PASS'; else; s = '**FAIL**'; end
end
