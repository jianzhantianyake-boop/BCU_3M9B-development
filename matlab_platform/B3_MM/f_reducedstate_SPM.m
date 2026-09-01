% =========================================================================
% 教学操作说明：稳定区域搜索中供 fsolve 调用的约化状态残差函数
% 使用方法：
%  按以下函数签名调用：function dfdt = f_reducedstate_SPM(x)
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
% x(2*nbus-ngen-1) = deltac(1: 2) | deltac_net(4:nbus) | V_net(10:9+nbus-ngen)
function dfdt = f_reducedstate_SPM(x)

postfault = evalin('base','postfault');
preset = evalin('base','preset');

delta2c=x(1);
delta3c=x(2);

Yfull = postfault.Yfull_mod;
ngen=size(preset.genno,1);
nbus=size(Yfull,1);
m=preset.m;
d=preset.d;
Pm=preset.Pmpu;
E=preset.Epu;
mT=sum(m,1);
Pe=zeros(1,ngen);
Pnet = zeros(nbus-ngen,1);
Qnet = zeros(nbus-ngen,1);
Transform = postfault.Transform;
G=real(Yfull);
B=imag(Yfull);

deltac = x(1:2);
delta1c= -m(2:ngen)'*deltac(1:2)/m(1);
deltacc = [delta1c delta2c delta3c];

deltac_net = [x(3) x(4) x(5) x(6) x(7) x(8)];
V_net=[x(9) x(10) x(11) x(12) x(13) x(14)];

%%  power calculation
% P calculation of gen
for i=1:ngen
    for j=1:ngen
        ddelta=deltacc(i)-deltacc(j);
        Pe(i)=Pe(i)+E(i)*E(j)*B(i,j)*sin(ddelta)+E(i)*E(j)*G(i,j)*cos(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=deltacc(i)-deltac_net(l);
        Pe(i)=Pe(i)+E(i)*V_net(l)*B(i,l+ngen)*sin(ddelta)+E(i)*V_net(l)*G(i,l+ngen)*cos(ddelta);
    end
end
% P calculation of Bus
for i=1:(nbus-ngen)
    for j=1:ngen
        ddelta=deltac_net(i)-deltacc(j);
        Pnet(i)=Pnet(i)+V_net(i)*E(j)*B(i+ngen,j)*sin(ddelta)+V_net(i)*E(j)*G(i+ngen,j)*cos(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=deltac_net(i)-deltac_net(l);
        Pnet(i)=Pnet(i)+V_net(i)*V_net(l)*B(i+ngen,l+ngen)*sin(ddelta)+V_net(i)*V_net(l)*G(i+ngen,l+ngen)*cos(ddelta);
    end
    for h=1:size(preset.Sload,1)
        if (preset.Sload(h,1)==Transform(i+ngen))
                   Pnet(i)=Pnet(i)+preset.Sload(h,2);
        end
    end
end
% Q calculation of Bus
for i=1:(nbus-ngen)
    for j=1:ngen
        ddelta=deltac_net(i)-deltacc(j);
        Qnet(i)=Qnet(i)-V_net(i)*E(j)*B(i+ngen,j)*cos(ddelta)+V_net(i)*E(j)*G(i+ngen,j)*sin(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=deltac_net(i)-deltac_net(l);
        Qnet(i)=Qnet(i)-V_net(i)*V_net(l)*B(i+ngen,l+ngen)*cos(ddelta)+V_net(i)*V_net(l)*G(i+ngen,l+ngen)*sin(ddelta);
    end
    for h=1:size(preset.Sload,1)
        if (preset.Sload(h,1)==Transform(i+ngen))
               Qnet(i)=Qnet(i)+preset.Sload(h,3);          
        end
    end
end

Pcoi=sum(Pm)-sum(Pe);

dfdt(1) =  (Pm(2)-Pe(2)-Pcoi/sum(m)*m(2))/m(2); 
dfdt(2) =  (Pm(3)-Pe(3)-Pcoi/sum(m)*m(3))/m(3);
% dfdt(1) =  (Pm(2)-Pe(2)-Pcoi/sum(m)*m(2))/d(2); 
% dfdt(2) =  (Pm(3)-Pe(3)-Pcoi/sum(m)*m(3))/d(3);

% power of load bus
dfdt(3) = -Pnet(1);
dfdt(4) = -Pnet(2);
dfdt(5) = -Pnet(3);
dfdt(6) = -Pnet(4);
dfdt(7) = -Pnet(5);
dfdt(8) = -Pnet(6);
% reactive power of load bus
dfdt(9) = -Qnet(1);
dfdt(10) = -Qnet(2);
dfdt(11) = -Qnet(3);
dfdt(12) = -Qnet(4);
dfdt(13) = -Qnet(5);
dfdt(14) = -Qnet(6);


dfdt = dfdt.';

end
