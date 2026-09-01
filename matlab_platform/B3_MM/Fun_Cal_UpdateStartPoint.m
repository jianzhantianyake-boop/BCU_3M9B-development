% =========================================================================
% 教学操作说明：轨迹迭代起点更新函数
% 使用方法：
%  按以下函数签名调用：function [thetac_update,flag_update]=Fun_Cal_UpdateStartPoint(thetac_lastpoint,preset,postfault)
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
%% Function: update start point along the ray of thetac_last and postfault SEP
function [thetac_update,flag_update]=Fun_Cal_UpdateStartPoint(thetac_lastpoint,preset,postfault)
%% Settings
    Y_red=postfault.Yred;
    ngen=size(Y_red,1);
    G_post=real(Y_red);
    B_post=imag(Y_red);
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    flag_update=0;
%% initialization
    Pe=zeros(ngen,1);
    len_ray=1e-3;
    dir_ray=zeros(1,ngen);  % vector from SEP_post to thetac_lastpoint
    thetac_SEP=postfault.SEP_delta;
    dir_ray=thetac_lastpoint-thetac_SEP;
    n_itermax=fix(2*norm(dir_ray)/len_ray);
    dir_ray=dir_ray/norm(dir_ray);
    n_iter=1;
    Ep_obsv=[0 0 0];

    Ep0=zeros(n_itermax,1);
    flag_position=0;    %1--init point is inside of boundary 2--outside
    for i=1:ngen
        for j=1:ngen
            ddelta=thetac_lastpoint(i)-thetac_lastpoint(j);
            Pe(i)=Pe(i)+E(i)*E(j)*B_post(i,j)*sin(ddelta)+E(i)*E(j)*G_post(i,j)*cos(ddelta);
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
    [Ep(1),Ep(2),Ep(3)]=Fun_Cal_PotentialEnergy(preset,postfault,thetac_SEP,thetac_SEP);
    Ep_obsv(2)=sum(Ep);
    thetac_act=thetac_SEP+len_ray*dir_ray;
%     thetac_act=thetac_SEP+len_ray*dir_ray;
    [Ep(1),Ep(2),Ep(3)]=Fun_Cal_PotentialEnergy(preset,postfault,thetac_SEP,thetac_act);
    Ep_obsv(3)=sum(Ep);
    
    while(n_iter~=-1)
        Ep_obsv(1)=Ep_obsv(2);
        Ep_obsv(2)=Ep_obsv(3);
        thetac_act=thetac_act+len_ray*dir_ray;
        [Ep(1),Ep(2),Ep(3)]=Fun_Cal_PotentialEnergy(preset,postfault,thetac_SEP,thetac_act);
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
        thetac_update=thetac_act-len_ray*dir_ray;
    else
%         error('No local maximum point found!');
        fprintf('No local maximum point found!');
        flag_update=0;
        thetac_update=thetac_lastpoint-len_ray*dir_ray;
    end
end
