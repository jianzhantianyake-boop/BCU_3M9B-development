% =========================================================================
% 教学操作说明：轨迹迭代起点更新函数
% 使用方法：
%  按以下函数签名调用：function [deltac_update,theta_update,voltage_update,flag_update]=Fun_Cal_UpdateStartPoint_SPM(deltac_lastpoint,theta_lastpoint,voltage_lastpoint,preset,postfault)
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
%% Function: update start point along the ray of deltac_last and postfault SEP
function [deltac_update,theta_update,voltage_update,flag_update]=Fun_Cal_UpdateStartPoint_SPM(deltac_lastpoint,theta_lastpoint,voltage_lastpoint,preset,postfault)
%% Settings
    Yfull=postfault.Yfull_mod;
    ngen=preset.ngen;
    nbus=preset.nbus;
    G_post=real(Yfull);
    B_post=imag(Yfull);
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    flag_update=0;
    Basevalue = evalin('base','Basevalue');
%% initialization
    Pe=zeros(ngen,1);
    len_ray=1e-3;


    deltac_SEP=postfault.SEP_delta;
    theta_SEP=postfault.net_delta;
    voltage_SEP=postfault.net_voltage;
    dir_ray=deltac_lastpoint-deltac_SEP;


    n_itermax=fix(2*norm(dir_ray)/len_ray);
    dir_ray=dir_ray/norm(dir_ray);

    n_iter=1;
    Ep_obsv=[0 0 0];

    Ep0=zeros(n_itermax,1);
    flag_position=0;    %1--init point is inside of boundary 2--outside

    for  i=1:ngen
        for j=1:ngen
            ddelta=deltac_lastpoint(i)-deltac_lastpoint(j);
            Pe(i)=Pe(i)+E(i)*E(j)*(G_post(i,j)*cos(ddelta)+B_post(i,j)*sin(ddelta));
        end
        for l=1:(nbus-ngen)
            ddelta=deltac_lastpoint(i)-theta_lastpoint(l);
            Pe(i)=Pe(i)+E(i)*voltage_lastpoint(l)*B_post(i,l+ngen)*sin(ddelta)+E(i)*voltage_lastpoint(l)*G_post(i,l+ngen)*cos(ddelta);
        end
    end
    
%% Check the direction (Use current point's vector)
    dotproduct=dir_ray'*(Pm-Pe-m/sum(m)*sum(Pm-Pe));
    if(dotproduct>0)
        flag_position=2;
%         len_ray=-len_ray;
    else
        flag_position=1;
    end
%% Search for the local maximum Ep along ray(and extension)
    [Ep(1),Ep(2),Ep(3),Ep(4),Ep(5)]=Fun_Cal_PotentialEnergy_SPM(preset,postfault,deltac_SEP,theta_SEP,voltage_SEP);
    %[Ep(1),Ep(2),Ep(3)]=Fun_Cal_PotentialEnergy(preset,postfault,deltac_SEP,deltac_SEP);
    Ep_obsv(2)=sum(Ep);
    deltac_act=deltac_SEP+len_ray*dir_ray;
    theta_est = theta_SEP;
    voltage_est = voltage_SEP;
    [theta_act,voltage_act,flag_iter,n_it,err] = Fun_AEiteration_SPM(theta_est,voltage_est,deltac_act,preset,Basevalue,"postfault",1e4,1e-12);
    if(flag_iter~=1)
        error('no AE solution');
    end
    [Ep(1),Ep(2),Ep(3),Ep(4),Ep(5)]=Fun_Cal_PotentialEnergy_SPM(preset,postfault,deltac_act,theta_act,voltage_act);
    %[Ep(1),Ep(2),Ep(3)]=Fun_Cal_PotentialEnergy(preset,postfault,deltac_SEP,deltac_act);
    Ep_obsv(3)=sum(Ep);
    
    while(n_iter~=-1)
        Ep_obsv(1)=Ep_obsv(2);
        Ep_obsv(2)=Ep_obsv(3);
        deltac_act=deltac_act+len_ray*dir_ray;
        theta_est = theta_act;
        voltage_est = voltage_act;
        [theta_act,voltage_act,flag_iter,n_it,err] = Fun_AEiteration_SPM(theta_est,voltage_est,deltac_act,preset,Basevalue,"postfault",1e4,1e-12);
        if(flag_iter~=1)
            error('no AE solution');
        end
        [Ep(1),Ep(2),Ep(3),Ep(4),Ep(5)]=Fun_Cal_PotentialEnergy_SPM(preset,postfault,deltac_act,theta_act,voltage_act);
        %[Ep(1),Ep(2),Ep(3)]=Fun_Cal_PotentialEnergy(preset,postfault,deltac_SEP,deltac_act);
        Ep_obsv(3)=sum(Ep);
        Ep0(n_iter)=Ep_obsv(3);
        n_iter=n_iter+1;
        if(Ep_obsv(2)>Ep_obsv(1)&&Ep_obsv(2)>Ep_obsv(3))
            n_iter=-1;
        elseif(n_iter>n_itermax)
            break;
        end
    end

    if(n_iter==-1)
        flag_update=1;
        deltac_update=deltac_act-len_ray*dir_ray;
        theta_est = theta_act;
        voltage_est = voltage_act;
        [theta_update,voltage_update,flag_iter,n_it,err] = Fun_AEiteration_SPM(theta_est,voltage_est,deltac_update,preset,Basevalue,"postfault",1e4,1e-10);
    else
        %error('No local maximum point found!');
        fprintf('No local maximum point found!');
        flag_update=0;
        deltac_update=deltac_lastpoint;
        theta_est=theta_lastpoint;
        voltage_est=voltage_lastpoint;
        [theta_update,voltage_update,flag_iter,n_it,err] = Fun_AEiteration_SPM(theta_est,voltage_est,deltac_update,preset,Basevalue,"postfault",1e4,1e-10);
    end
end
