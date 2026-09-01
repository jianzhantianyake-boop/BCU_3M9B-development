function [ok, msgs] = bcu_validate_config(cfg)
% =========================================================================
% bcu_validate_config —— 运行前配置校验（fail-fast）
%
% 在真正跑仿真前把「维度不匹配、支路不存在、步长非法」等错误挡在门口，
% 避免跑到一半才崩，并给出可读的错误定位。返回：
%   ok   —— 逻辑值，是否全部通过
%   msgs —— cellstr，逐条校验结果（含 [OK]/[ERR]/[WARN] 前缀）
% =========================================================================
msgs = {};
nerr = 0;

    function add(level, txt)
        msgs{end+1} = sprintf('[%s] %s', level, txt); %#ok<AGROW>
        if strcmp(level,'ERR'), nerr = nerr + 1; end
    end

% --- 1. mode 合法性 ---
validModes = ["reduced_cct","reduced_numerical","reduced_region", ...
              "spm_cct","spm_numerical","spm_region", ...
              "two_machine_region_3d","two_machine_region_3d_gfl"];
if any(cfg.mode == validModes)
    add('OK', sprintf('mode = %s', cfg.mode));
else
    add('ERR', sprintf('未知 mode: %s', cfg.mode));
end

% --- 2. case 文件存在 ---
if exist(cfg.CaseName, 'file') == 2
    add('OK', sprintf('CaseName = %s (found)', cfg.CaseName));
else
    add('ERR', sprintf('找不到 case 文件: %s（检查 MATPOWER data 路径是否已加入 path）', cfg.CaseName));
end

% --- 3. 发电机向量维度一致 ---
ng = numel(cfg.m);
lens = [numel(cfg.m), numel(cfg.damping_ratio), numel(cfg.Pm), numel(cfg.xd1), numel(cfg.E)];
if all(lens == ng)
    add('OK', sprintf('发电机向量维度一致 (ngen = %d)', ng));
else
    add('ERR', sprintf('发电机向量长度不一致: m/damp/Pm/xd1/E = [%s]', num2str(lens)));
end

% --- 4. 与 case 实际发电机数匹配（能加载 case 时）---
if exist(cfg.CaseName,'file') == 2
    try
        mpc = feval(cfg.CaseName);
        ng_case = size(mpc.gen, 1);
        if ng_case == ng
            add('OK', sprintf('发电机数与 %s 匹配 (%d)', cfg.CaseName, ng));
        else
            add('ERR', sprintf('%s 有 %d 台机，但参数给了 %d 台', cfg.CaseName, ng_case, ng));
        end
        % --- 5. 故障支路存在性 ---
        fl = cfg.faultline(:);
        br = mpc.branch(:, 1:2);
        hit = any((br(:,1)==fl(1) & br(:,2)==fl(2)) | (br(:,1)==fl(2) & br(:,2)==fl(1)));
        if hit
            add('OK', sprintf('故障支路 [%d-%d] 存在于 case', fl(1), fl(2)));
        else
            add('ERR', sprintf('故障支路 [%d-%d] 在 %s 中不存在', fl(1), fl(2), cfg.CaseName));
        end
    catch ME
        add('WARN', sprintf('无法加载 case 做深度校验: %s', ME.message));
    end
end

% --- 6. 故障位置 ---
if any(cfg.faultposition == [0 1])
    add('OK', sprintf('faultposition = %d', cfg.faultposition));
else
    add('ERR', sprintf('faultposition 必须为 0 或 1，当前 %g', cfg.faultposition));
end

% --- 7. 求解器 / 数值参数 ---
if any(cfg.EquCal == [1 2]); add('OK','EquCal 合法'); else; add('ERR','EquCal 必须为 1 或 2'); end
if cfg.PathEnergyCal >= -1 && cfg.PathEnergyCal == round(cfg.PathEnergyCal)
    add('OK','PathEnergyCal 合法');
else
    add('ERR','PathEnergyCal 必须为 -1 或非负整数');
end
if cfg.Tunit > 0 && cfg.Tfault > 0 && cfg.Tunit < cfg.Tfault
    add('OK', sprintf('Tfault=%.3g s, Tunit=%.1e s (%d 步)', cfg.Tfault, cfg.Tunit, round(cfg.Tfault/cfg.Tunit)));
else
    add('ERR', 'Tfault/Tunit 非法：需 Tunit>0 且 0<Tunit<Tfault');
end
if cfg.f_base > 0; add('OK', sprintf('f_base = %g Hz', cfg.f_base)); else; add('ERR','f_base 必须为正'); end

% --- 8. 物理合理性（软警告）---
sumPm = sum(cfg.Pm);
if sumPm <= 0
    add('WARN', sprintf('ΣPm = %.3f <= 0，潮流/平衡点可能异常', sumPm));
end

ok = (nerr == 0);
end
