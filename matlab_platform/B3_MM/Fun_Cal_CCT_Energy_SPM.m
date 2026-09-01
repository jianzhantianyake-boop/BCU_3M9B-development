% =========================================================================
% 教学操作说明：临界切除时间或临界能量相关量的计算函数
% 使用方法：
%  按以下函数签名调用：function [CCT,Exit_deltac,Exit_omegac,Exit_theta,Exit_omega,Exit_voltage,flag_CCT]=Fun_Cal_CCT_Energy_SPM(E_critical,fault,postfault,preset)
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
%% function: calculate CCT along fault-on trajectory where Ek+Ep=Ecritical
function [CCT,Exit_deltac,Exit_omegac,Exit_theta,Exit_omega,Exit_voltage,flag_CCT]=Fun_Cal_CCT_Energy_SPM(E_critical,fault,postfault,preset)
    cycle=size(fault.traj.omega,1);
    Ek=zeros(cycle,1);
    Ep=zeros(cycle,1);
    Esum=zeros(cycle,1);
    ngen=size(postfault.Yred,1);
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    flag_CCT=0;
    CCT=0;
    for tm=1:cycle
        for i=1:ngen
        Ek(tm)=Ek(tm)+0.5*m(i)*fault.traj.omegac(tm,i)^2;
        end
        [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3),Ep_tmp(4),Ep_tmp(5)]=Fun_Cal_PotentialEnergy_SPM(preset,postfault,fault.traj.deltac(tm,:)',fault.traj.theta(tm,:)',fault.traj.voltage(tm,:)');
        Ep(tm)=sum(Ep_tmp);
        clear Ep_tmp
        Esum(tm)=Ep(tm)+Ek(tm);
        if(tm>1)
        if(Esum(tm-1)<E_critical&&Esum(tm)>E_critical&&flag_CCT==0)
            CCT=(tm-1)*fault.traj.Tunit;
            Exit_deltac=fault.traj.deltac(tm-1,:);
            Exit_omegac=fault.traj.omegac(tm-1,:);
            Exit_theta=fault.traj.theta(tm-1,:);
            Exit_omega=fault.traj.omega(tm-1,:);
            Exit_voltage = fault.traj.voltage(tm-1,:);
            flag_CCT=1;
        end
        end
    end
end
