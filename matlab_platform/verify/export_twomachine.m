% =========================================================================
% export_twomachine.m -- 导出 two_machine_region_3d 与 _gfl 的平衡点参考
% 注: Statable_Region_3D*.m 顶部有 `clear;`, 会清普通变量; 故路径用 global 保存,
%     run 后重新声明恢复; 提取用 local function(不受 clear 影响).
% =========================================================================
global VDIR VROOT
VDIR  = fileparts(mfilename('fullpath'));
VROOT = fileparts(VDIR);
addpath(VROOT); setup_bcu_paths(VROOT);
set(0, 'DefaultFigureVisible', 'off');

cd(VROOT);
run(fullfile('B3_MM', 'Statable_Region_3D.m'));
global VDIR VROOT
save_eps(ep_set, fullfile(VDIR, 'baseline_twomachine.mat'), 'TM');

cd(VROOT);
run(fullfile('B3_MM', 'Statable_Region_3D_GFL.m'));
global VDIR VROOT
save_eps(ep_set, fullfile(VDIR, 'baseline_twomachine_gfl.mat'), 'TMGFL');

function save_eps(ep_set, outfile, tag)
    nep = numel(ep_set);
    xeps = cell(nep, 1); flags = zeros(nep, 1);
    for i = 1:nep
        xeps{i} = double(ep_set(i).xep(:)).';
        flags(i) = ep_set(i).flag;
    end
    dims = cellfun(@numel, xeps);
    save(outfile, 'xeps', 'flags', '-v7');
    fprintf('EXPORT_%s_OK nep=%d dims=[%s]\n', tag, nep, num2str(dims(:).'));
end
