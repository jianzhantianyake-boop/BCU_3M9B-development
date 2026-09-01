% =========================================================================
% 教学操作说明：稳定区域搜索中供 fsolve 调用的约化状态残差函数
% 使用方法：
%  按以下函数签名调用：function dfdt = f_reducedstate(deltac)
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
function dfdt = f_reducedstate(deltac)

postfault = evalin('base','postfault');
preset = evalin('base','preset');

delta2c=deltac(1);
delta3c=deltac(2);

Y_red=postfault.Yred;
G=real(Y_red);    
B=imag(Y_red);
ngen=size(Y_red,1);
m=preset.m;
d=preset.d;
Pm=preset.Pmpu;
E=preset.Epu;
mT=sum(m,1);
Pe=zeros(1,ngen);

delta1c= -m(2:ngen)'*deltac(1:2)/m(1);
deltacc = [delta1c delta2c delta3c];

% Pe calculation
for i=1:ngen
    for j=1:ngen
        ddelta=deltacc(i)-deltacc(j);
        Pe(i)=Pe(i)+E(i)*E(j)*(G(i,j)*cos(ddelta)+B(i,j)*sin(ddelta));
    end
end
Pcoi=sum(Pm)-sum(Pe);

dfdt(1) =  (Pm(2)-Pe(2)-Pcoi/sum(m)*m(2)); 
dfdt(2) =  (Pm(3)-Pe(3)-Pcoi/sum(m)*m(3));


dfdt = dfdt.';

end
