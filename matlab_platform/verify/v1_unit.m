% =========================================================================
% v1_unit.m —— 第 1 层：底层函数的解析解/守恒律单元测试
%
% 依赖 v0_baseline.m 先生成的 baseline_reduced.mat（提供 preset/pre-/postfault
% /fault/Basevalue 等真实参数）。所有断言用 [[UNIT]] 前缀打印，英文输出避免
% GBK/UTF-8 乱码。
%
% 覆盖：
%   U1 摆动方程右端在 SEP 处应为 0（平衡点定义）
%   U2 势能函数起点=终点时三项应为 0（路径无关基准）
%   U3 无损网络(G=0)时路径相关势能项 Ep3 应为 0
%   U4 势能路径积分自洽：Ray近似 vs 多段梯形积分应收敛
%   U5 RK4 步长收敛：步长减半，退出点/CCT 应几乎不变（检验积分器阶数）
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisDir);
addpath(projectRoot); setup_bcu_paths(projectRoot);
load(fullfile(thisDir,'baseline_reduced.mat'), ...
     'preset','prefault','postfault','fault','Basevalue');

npass = 0; nfail = 0;
fprintf('\n========== [[UNIT]] 第 1 层单元测试 ==========\n');

% ---- U1: f_reducedstate 在 postfault SEP 处右端应为 0 ----
% f_reducedstate 通过 evalin('base',...) 取 postfault/preset，需在 base。
sep = postfault.SEP_delta;
r = f_reducedstate([sep(2); sep(3)]);
[npass,nfail] = chk('U1 swing RHS at SEP ~ 0', norm(r) < 1e-6, ...
    sprintf('|f(SEP)|=%.3e', norm(r)), npass, nfail);

% ---- U2: 势能起点=终点 => 三项全 0 ----
[e1,e2,e3] = Fun_Cal_PotentialEnergy(preset,postfault, ...
    postfault.SEP_delta, postfault.SEP_delta);
[npass,nfail] = chk('U2 PE(SEP,SEP) == 0', max(abs([e1 e2 e3])) < 1e-12, ...
    sprintf('[Ep1 Ep2 Ep3]=[%.2e %.2e %.2e]', e1,e2,e3), npass, nfail);

% ---- U3: 无损网络 (G=0) => 路径项 Ep3 = 0 ----
pf_lossless = postfault;
pf_lossless.Yred = 1i*imag(postfault.Yred);   % 去掉电导，仅留电纳
endpt = postfault.CUEP_delta;                 % 任取非平凡终点
[~,~,e3_lossless] = Fun_Cal_PotentialEnergy(preset, pf_lossless, ...
    postfault.SEP_delta, endpt);
[npass,nfail] = chk('U3 lossless => Ep3 == 0', abs(e3_lossless) < 1e-12, ...
    sprintf('Ep3_lossless=%.3e', e3_lossless), npass, nfail);

% ---- U4: 路径积分自洽 —— Ray近似(0) vs 多段梯形(大N) 应接近 ----
endpt = postfault.CUEP_delta;
pr0 = preset; pr0.PathEnergyCal = 0;      % Ray approximation
[~,~,e3_ray]  = Fun_Cal_PotentialEnergy(pr0, postfault, postfault.SEP_delta, endpt);
prN = preset; prN.PathEnergyCal = 200;    % 200 段梯形积分（细分路径）
[~,~,e3_trap] = Fun_Cal_PotentialEnergy(prN, postfault, postfault.SEP_delta, endpt);
rel = abs(e3_ray - e3_trap) / max(abs(e3_trap), 1e-12);
[npass,nfail] = chk('U4 path-integral consistency (Ray vs trap)', rel < 0.05, ...
    sprintf('Ep3_ray=%.4e Ep3_trap=%.4e rel=%.2f%%', e3_ray, e3_trap, 100*rel), ...
    npass, nfail);

% ---- U5: RK4 步长收敛 —— 步长减半，退出点应几乎不变 ----
delta0 = prefault.SEP_delta;
omega0 = prefault.SEP_omegapu * Basevalue.omegab;
Tf = 0.35;                                 % 覆盖到 CCT(~0.24) 之后即可
[~,~,thc_a,~,tm_a] = Fun_Cal_Exitpoint(Tf,1e-4,   fault.Yred,postfault.Yred,preset,delta0,omega0,Basevalue.omegab);
[~,~,thc_b,~,tm_b] = Fun_Cal_Exitpoint(Tf,0.5e-4, fault.Yred,postfault.Yred,preset,delta0,omega0,Basevalue.omegab);
exit_a = thc_a(tm_a,:);
exit_b = thc_b(tm_b,:);
t_exit_a = (tm_a-1)*1e-4;  t_exit_b = (tm_b-1)*0.5e-4;
d_exit = norm(exit_a - exit_b);
[npass,nfail] = chk('U5 RK4 step-halving convergence', d_exit < 1e-3, ...
    sprintf('t_exit: %.4f vs %.4f s | |dExit|=%.3e rad', t_exit_a, t_exit_b, d_exit), ...
    npass, nfail);

fprintf('---------- [[UNIT]] 通过 %d / 失败 %d ----------\n', npass, nfail);

% ------------------------- helpers --------------------------------
function [np,nf] = chk(name, cond, detail, np, nf)
    if cond
        fprintf('[[UNIT]] PASS  %-42s | %s\n', name, detail); np = np+1;
    else
        fprintf('[[UNIT]] FAIL  %-42s | %s\n', name, detail); nf = nf+1;
    end
end
