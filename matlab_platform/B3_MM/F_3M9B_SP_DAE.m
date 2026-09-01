% =========================================================================
% 教学操作说明：三机九母线模型的状态导数计算函数
% 使用方法：
%  按以下函数签名调用：function dfdt = F_3M9B_SP_DAE(x,system)
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
% x(2*nbus-ngen) = deltac(1: 3) | deltac_net(4:nbus) |
% V_net(10:9+nbus-ngen)|omega(16:18)
function dfdt = F_3M9B_SP_DAE(x,system)
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

deltac = [x(1) x(2) x(3)];
omega1 = x(16);
omega2 = x(17);
omega3 = x(18);
omega_coi = [omega1 omega2 omega3]*m/mT;


%system = evalin('base','system');
switch system
    case "prefault"
        Yfull = prefault.Yfull_mod;
        Transform = prefault.Transform;
    case "fault1"
        Yfull = fault.Yfull_mod;
        Transform = fault.Transform;
    case "fault2"
        Yfull = fault.Yfull_mod2;
        Transform = fault.Transform2;
    case "postfault"
        Yfull = postfault.Yfull_mod;
        Transform = postfault.Transform;
end

if system == "fault1"
    deltac_net = [x(4) x(5) x(6) x(7) x(8)];
    V_net=[x(10) x(11) x(12) x(13) x(14)];
else
    deltac_net = [x(4) x(5) x(6) x(7) x(8) x(9)];
    V_net=[x(10) x(11) x(12) x(13) x(14) x(15)];
end

G=real(Yfull);
B=imag(Yfull);
ngen=size(preset.genno,1);
nbus=size(Yfull,1);
Pe=zeros(1,ngen);
Pnet = zeros(nbus-ngen,1);
Qnet = zeros(nbus-ngen,1);

%%  power calculation
% P calculation of gen
for i=1:ngen
    for j=1:ngen
        ddelta=deltac(i)-deltac(j);
        Pe(i)=Pe(i)+E(i)*E(j)*B(i,j)*sin(ddelta)+E(i)*E(j)*G(i,j)*cos(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=deltac(i)-deltac_net(l);
        Pe(i)=Pe(i)+E(i)*V_net(l)*B(i,l+ngen)*sin(ddelta)+E(i)*V_net(l)*G(i,l+ngen)*cos(ddelta);
    end
end
% P calculation of Bus
for i=1:(nbus-ngen)
    for j=1:ngen
        ddelta=deltac_net(i)-deltac(j);
        Pnet(i)=Pnet(i)+V_net(i)*E(j)*B(i+ngen,j)*sin(ddelta)+V_net(i)*E(j)*G(i+ngen,j)*cos(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=deltac_net(i)-deltac_net(l);
        Pnet(i)=Pnet(i)+V_net(i)*V_net(l)*B(i+ngen,l+ngen)*sin(ddelta)+V_net(i)*V_net(l)*G(i+ngen,l+ngen)*cos(ddelta);
    end
    for h=1:size(preset.Sload,1)
        if (preset.Sload(h,1)==Transform(i+ngen))
               if (system == "fault1")||(system == "fault2")
                   if Transform(i+ngen)~=fault.faultbus
                        %Pnet(i)=Pnet(i)+preset.Sload(h,2);% during fault
                        %pure impedance
                   end
               else
                   Pnet(i)=Pnet(i)+preset.Sload(h,2);
               end
        end
    end
end


% Q calculation of Bus
for i=1:(nbus-ngen)
    for j=1:ngen
        ddelta=deltac_net(i)-deltac(j);
        Qnet(i)=Qnet(i)-V_net(i)*E(j)*B(i+ngen,j)*cos(ddelta)+V_net(i)*E(j)*G(i+ngen,j)*sin(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=deltac_net(i)-deltac_net(l);
        Qnet(i)=Qnet(i)-V_net(i)*V_net(l)*B(i+ngen,l+ngen)*cos(ddelta)+V_net(i)*V_net(l)*G(i+ngen,l+ngen)*sin(ddelta);
    end
    for h=1:size(preset.Sload,1)
        if (preset.Sload(h,1)==Transform(i+ngen))
           if (system == "fault1")||(system == "fault2")
               if Transform(i+ngen)~=fault.faultbus
                    %Qnet(i)=Qnet(i)+preset.Sload(h,3);% during fault
                        %pure impedance
               end
           else
               Qnet(i)=Qnet(i)+preset.Sload(h,3);
           end
               
        end
    end
end



%% differential equation
% deltac
dfdt(1) = omega1 - omega_coi;% deltac1
dfdt(2) = omega2 - omega_coi;% deltac2
dfdt(3) = omega3 - omega_coi;% deltac3

%omega
dfdt(16) = (Pm(1)-Pe(1)-d(1)*(omega1-omegab))/m(1);% omega1
dfdt(17) = (Pm(2)-Pe(2)-d(2)*(omega2-omegab))/m(2);% omega2
dfdt(18) = (Pm(3)-Pe(3)-d(3)*(omega3-omegab))/m(3);% omega3

%% agebric equation
if system == "fault1"
% power of load bus
dfdt(4) = -Pnet(1);
dfdt(5) = -Pnet(2);
dfdt(6) = -Pnet(3);
dfdt(7) = -Pnet(4);
dfdt(8) = -Pnet(5);
dfdt(9) = 0;
% reactive power of load bus
dfdt(10) = -Qnet(1);
dfdt(11) = -Qnet(2);
dfdt(12) = -Qnet(3);
dfdt(13) = -Qnet(4);
dfdt(14) = -Qnet(5);
dfdt(15) = 0;
else
% power of load bus
dfdt(4) = -Pnet(1);
dfdt(5) = -Pnet(2);
dfdt(6) = -Pnet(3);
dfdt(7) = -Pnet(4);
dfdt(8) = -Pnet(5);
dfdt(9) = -Pnet(6);
% reactive power of load bus
dfdt(10) = -Qnet(1);
dfdt(11) = -Qnet(2);
dfdt(12) = -Qnet(3);
dfdt(13) = -Qnet(4);
dfdt(14) = -Qnet(5);
dfdt(15) = -Qnet(6);
end

dfdt = dfdt.';

end
