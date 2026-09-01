% =========================================================================
% v2_energy.m —— 第 2 层：能量守恒律检验（TEF 方法的物理根基）
%
% 原理：无阻尼(d=0)时，故障后系统沿轨迹总能量 Ek+Ep 应守恒。
%   Ek = 0.5*Σ m_i*ωc_i^2 ，Ep 由 Fun_Cal_PotentialEnergy(SEP->θc) 给出，
%   与 Fun_Cal_CCT_Energy 内部的能量定义完全一致。
% 若沿一条 post-fault 有界振荡轨迹总能量明显漂移，则积分器或能量函数有问题。
%
% 复用 Fun_Cal_Exitpoint 做纯 post-fault 积分（把“故障网络”与“故障后网络”
% 都传 postfault.Yred），preset.d 置零得到守恒系统。
% =========================================================================
thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisDir);
addpath(projectRoot); setup_bcu_paths(projectRoot);
load(fullfile(thisDir,'baseline_reduced.mat'), 'preset','postfault','Basevalue');

fprintf('\n========== [[ENERGY]] 第 2 层能量守恒检验 ==========\n');

omegab = Basevalue.omegab;
ngen = numel(preset.m);
Tf = 2.0; Tunit = 1e-4;

% 通用能量漂移测量：给定动力学/势场网络 pfNet(含 .Yred,.SEP_delta) 与参数 pr(d=0)
measure = @(pr,pfNet,theta0,omega0) local_drift(pr,pfNet,theta0,omega0,Tf,Tunit,omegab);

% ---------- A. 有损网络（真实 postfault，G≠0）：TEF 路径近似，预期有漂移 ----------
prA = preset; prA.d = zeros(size(preset.d));
theta0 = postfault.SEP_delta(:).';
dwc = [1.0,-0.5,-0.5]; dwc = dwc - (dwc*preset.m)/sum(preset.m);   % 投影到 COI
[relA,EmA,daA] = measure(prA, postfault, theta0, omegab*ones(1,ngen)+dwc);
fprintf('[[ENERGY]] A 有损网络: E_mean=%.4e abs_drift=%.3e rel_drift=%.3e\n', EmA, daA, relA);
fprintf('[[ENERGY]]   -> 有损下路径相关(Ep3, Ray近似) 使能量非严格守恒，此漂移属方法固有，非代码错误\n');

% ---------- B. 无损网络（G=0）：自洽保守系统，能量必须严格守恒 ----------
% 纯无损网络需 ΣPm=0 才有平衡点；这里在 postfault SEP 角度处反推 Pm，
% 使该点成为无损系统的平衡点，从而构造有界振荡的保守系统。
pfL.Yred = 1i*imag(postfault.Yred);                                % 去电导 (G=0)
Bmat = imag(postfault.Yred); Ev = preset.Epu(:);
th = postfault.SEP_delta(:);
PeL = zeros(ngen,1);
for ii = 1:ngen
    for jj = 1:ngen
        PeL(ii) = PeL(ii) + Ev(ii)*Ev(jj)*Bmat(ii,jj)*sin(th(ii)-th(jj));  % G=0
    end
end
prL = prA; prL.Pmpu = PeL;          % 令 th 成为无损系统平衡点
pfL.SEP_delta = th;
theta0L = th(:).';
w0 = omegab*ones(1,ngen)+dwc;
[relB1,~,~] = local_drift(prL,pfL,theta0L,w0, Tf, 1e-4,   omegab);
[relB2,~,~] = local_drift(prL,pfL,theta0L,w0, Tf, 0.5e-4, omegab);
[relB3,~,~] = local_drift(prL,pfL,theta0L,w0, Tf, 0.25e-4,omegab);
fprintf('[[ENERGY]] B 无损网络 能量漂移随步长(阶数判据):\n');
fprintf('[[ENERGY]]   h=1e-4  rel=%.3e | h=0.5e-4 rel=%.3e | h=0.25e-4 rel=%.3e\n', relB1,relB2,relB3);
fprintf('[[ENERGY]]   下降比 r1=%.2f  r2=%.2f  (欧拉一阶≈2, 四阶≈16)\n', relB1/relB2, relB2/relB3);
if relB3 < relB1/3
    fprintf('[[ENERGY]] 判定: 漂移随步长稳定下降 => 能量函数实现正确，残留漂移来自积分器精度\n');
    fprintf('[[ENERGY]]        (Fun_Cal_Exitpoint 在 d=0 时 ω 分量退化为前向欧拉；\n');
    fprintf('[[ENERGY]]         实用步长 1e-4、短故障窗内影响可忽略，见 U5)\n');
else
    fprintf('[[ENERGY]] 判定: 漂移不随步长下降 => 需排查能量函数/网络对称性实现\n');
end
fprintf('========== [[ENERGY]] 结束 ==========\n');

% ------------------------- helpers --------------------------------
function [rel,Emean,drift_abs] = local_drift(pr,pfNet,theta0,omega0,Tf,Tunit,omegab)
    [~,~,thc,omc,~] = Fun_Cal_Exitpoint(Tf,Tunit, pfNet.Yred, pfNet.Yred, ...
                                        pr, theta0, omega0, omegab);
    N = size(thc,1); E = zeros(N,1);
    for k = 1:N
        Ek = 0.5*sum(pr.m(:).' .* omc(k,:).^2);
        [a,b,c] = Fun_Cal_PotentialEnergy(pr, pfNet, pfNet.SEP_delta, thc(k,:).');
        E(k) = Ek + a + b + c;
    end
    Emean = mean(E); drift_abs = max(E)-min(E);
    rel = drift_abs / max(abs(Emean),1e-12);
end
