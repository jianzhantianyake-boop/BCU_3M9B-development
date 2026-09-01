% =========================================================================
% 教学操作说明：三机九母线模型的状态导数计算函数
% 使用方法：
%  按以下函数签名调用：function dfdt = F_3M9B_MR_ODE(deltacomega)
% 参数：
%  以签名中的输入变量为准。状态向量、导纳矩阵和参数结构体的维度及排序必须与初始化结果一致。
% 返回 / 工作区结果：
%  以签名左侧的输出变量为准；结果通常供上层 CCT、轨迹、能量或稳定区域脚本继续使用。
% 步骤：
%  1. 接收上层脚本提供的状态与参数。 2. 按原始方程或迭代规则完成本步骤。 3. 返回残差、状态、能量、矩阵或收敛标志。
% 单位：
%   角度通常为 rad，角速度为 rad/s，时间为 s，功率、电压和导纳通常为 pu；
%   个别中间变量为无量纲标志、迭代次数或矩阵索引，具体以变量定义为准。
% 前置条件：
%  仅在其上层流程已完成初始化后调用；不要仅凭单次函数输出推断整个系统的稳定性。
% 研究与验证边界：
%   本次只更新教学注释，不改变原始方程、参数、判据、求解器或绘图逻辑；
%   MATLAB 原生运行、收敛性和物理结论仍须在目标 MATLAB 环境中实际核验。
% =========================================================================
function dfdt = F_3M9B_MR_ODE(deltacomega)
%% parameters preprocess  
preset = evalin('base','preset');
prefault = evalin('base','prefault');
fault = evalin('base','fault');
postfault = evalin('base','postfault');
Basevalue = evalin('base','Basevalue');
Pm=preset.Pmpu;
E=preset.Epu;
m=preset.m;
d=preset.d;
mT=sum(m,1);
omegab = Basevalue.omegab;

deltac = [deltacomega(1) deltacomega(2) deltacomega(3)];
omega1 = deltacomega(4);
omega2 = deltacomega(5);
omega3 = deltacomega(6);
omega_coi = [omega1 omega2 omega3]*m/mT;


system = evalin('base','system');
switch system
    case "prefault"
        Yred = prefault.Yred;
    case "fault"
        Yred = fault.Yred;
    case "postfault"
        Yred = postfault.Yred;
end

G=real(Yred);
B=imag(Yred);
ngen=size(Yred,1);
Pe=zeros(1,ngen);

%% power calculation
for i=1:ngen
    for j=1:ngen
        delta=deltac(i)-deltac(j);
        Pe(i)=Pe(i)+E(i)*E(j)*(G(i,j)*cos(delta)+B(i,j)*sin(delta));
    end
end


%% differential equation
% deltac
dfdt(1) = omega1 - omega_coi;% deltac1
dfdt(2) = omega2 - omega_coi;% deltac2
dfdt(3) = omega3 - omega_coi;% deltac3

%omega
dfdt(4) = (Pm(1)-Pe(1)-d(1)*(omega1-omegab))/m(1);% omega1
dfdt(5) = (Pm(2)-Pe(2)-d(2)*(omega2-omegab))/m(2);% omega2
dfdt(6) = (Pm(3)-Pe(3)-d(3)*(omega3-omegab))/m(3);% omega3


dfdt = dfdt.';

end
