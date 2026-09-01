% =========================================================================
% 教学操作说明：稳定平衡点误差与收敛状态检查函数
% 使用方法：
%  按以下函数签名调用：function [Perr,flag_err]=Fun_SEPcheck(State,preset,delta_SEP,omega_SEP)
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
%% function: check the err of SEP result
%% input: omega_SEP (in SRF frame, rad/s)
function [Perr,flag_err]=Fun_SEPcheck(State,preset,delta_SEP,omega_SEP)
    Yred=State.Yred;
    ngen=size(Yred,1);
    B=imag(Yred);
    G=real(Yred);
    E=preset.Epu;
    d=preset.d;
    m=preset.m;
    Pm=preset.Pmpu;
    flag_err=0;
    Pe_tmp=zeros(ngen,1);
    for i=1:ngen
        for j=1:ngen
            ddelta=delta_SEP(i)-delta_SEP(j);
            Pe_tmp(i)=Pe_tmp(i)+E(i)*E(j)*(B(i,j)*sin(ddelta)+G(i,j)*cos(ddelta));
        end
    end
    Pcoi_tmp=sum(Pm)-sum(Pe_tmp)-sum(d)*omega_SEP;
    Err_PCUEP=Pm-Pe_tmp-Pcoi_tmp/sum(m)*m-d*omega_SEP;
    Perr=Err_PCUEP;
    if(norm(Perr)>1e-2)
        flag_err=1;
    end
end
