% =========================================================================
% 教学操作说明：临界切除时间或临界能量相关量的计算函数
% 使用方法：
%  按以下函数签名调用：function [CCT,Exit_thetac,Exit_omegac,Exit_theta,Exit_omega,flag_CCT,Traj_Stb,Traj_Unstb]=Fun_Cal_CCT_Real(fault,postfault,preset,Basevalue,CCT_ref)
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
%% function: calculate real CCT along fault-on trajectory
%% unstable criterion: 
function [CCT,Exit_thetac,Exit_omegac,Exit_theta,Exit_omega,flag_CCT,Traj_Stb,Traj_Unstb]=Fun_Cal_CCT_Real(fault,postfault,preset,Basevalue,CCT_ref)
    cycle=size(fault.traj.omega,1);
    ngen=size(postfault.Yred,1);
    omegab=Basevalue.omegab;
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    flag_CCT=0;
    Tunitmin=fault.traj.Tunit;
    Tfaultlen=fault.traj.Tlength;
    thetac_fault=fault.traj.thetac;
    omega_fault=fault.traj.omega;
    
%% Change unit settings
    Changeunit.n_stepchange=5; % numbers of units
    Changeunit.candidate=[1e3;1e2;10;5;1];
    if(Tfaultlen<=2*Tunitmin*Changeunit.candidate(1))
        error('Unsuitable Tunit selected, choose smaller T_unit or regenerate fault traj with longer time');
    end
%% Stable & Unstable settings
    Stb.n_current=fix(CCT_ref/Tunitmin);    % start from reference CCT
    Stb.n_stbmax=1;
    Stb.n_unstbmin=cycle;
%% Postfault traj iteration settings
    Tunit=1e-4;
    Tpostmax=200;    % maximum time for postfault iteration
%% Prejudge the CCT_ref
    delta0=thetac_fault(Stb.n_current,:);
    omega0=omega_fault(Stb.n_current,:);
    [theta_RK4,omega_RK4,thetac_RK4,omegacoi_RK4,Pe,cycle,flag_unstb]=Fun_TrajIter_StableCheck_SRF(Tpostmax,Tunit,postfault,preset,delta0,omega0,omegab);
    if(flag_unstb==0)
        Stb.n_stbmax=Stb.n_current;
        Stb.n_current=(fix(Stb.n_stbmax/Changeunit.candidate(1))+1)*Changeunit.candidate(1);
        Traj_Stb.thetac=fault.traj.thetac;
        Traj_Stb.omegac=fault.traj.omegac;
        Traj_Stb.theta=fault.traj.theta;
        Traj_Stb.omega=fault.traj.omega;
        
    else
        Stb.n_unstbmin=Stb.n_current;
        Stb.n_current=(fix(Stb.n_unstbmin/Changeunit.candidate(1)))*Changeunit.candidate(1);
    end
%% Change unit
    for n_step=1:Changeunit.n_stepchange
        Changeunit.Currentunit=Changeunit.candidate(n_step);
        flag_changeunit=0;
    %% Current unit scan
    while(flag_changeunit==0)
        delta0=thetac_fault(Stb.n_current,:);
        omega0=omega_fault(Stb.n_current,:);
        [theta_RK4,omega_RK4,thetac_RK4,omegacoi_RK4,Pe,cycle,flag_unstb]=Fun_TrajIter_StableCheck_SRF(Tpostmax,Tunit,postfault,preset,delta0,omega0,omegab);
        
        if(flag_unstb==0)
            if(Stb.n_stbmax<Stb.n_current)
                Stb.n_stbmax=Stb.n_current;
                Traj_Stb.thetac=thetac_RK4(1:cycle,:);
                Traj_Stb.omegac=omega_RK4-omegacoi_RK4*ones(1,ngen);
                Traj_Stb.theta=theta_RK4(1:cycle,:);
                Traj_Stb.omega=omega_RK4(1:cycle,:);
            end
        else
            if(Stb.n_unstbmin>Stb.n_current)    % for initialization
                Stb.n_unstbmin=Stb.n_current;                    
                Traj_Unstb.thetac=thetac_RK4(1:cycle,:);
                Traj_Unstb.omegac=omega_RK4(1:cycle,:)-omegacoi_RK4(1:cycle)*ones(1,ngen);
                Traj_Unstb.theta=theta_RK4(1:cycle,:);
                Traj_Unstb.omega=omega_RK4(1:cycle,:);
            end
        end
        if(Stb.n_unstbmin-Stb.n_stbmax<=Changeunit.Currentunit)
            flag_changeunit=1;
            if(n_step<Changeunit.n_stepchange)
                Changeunit.Nextunit=Changeunit.candidate(n_step+1);
                Stb.n_current=(fix(Stb.n_stbmax/Changeunit.Nextunit)+1)*Changeunit.Nextunit;
            else
                CCT=Stb.n_stbmax*Changeunit.candidate(n_step)*Tunitmin;
                flag_CCT=1;
            end
        end
        if(Stb.n_unstbmin<=Stb.n_stbmax)
            error('Stb.n_unstbmin<=Stb.n_stbmax');
        end
        if(flag_changeunit==0)
            if(flag_unstb==0)
                Stb.n_current=Stb.n_current+Changeunit.Currentunit;
            else
                Stb.n_current=Stb.n_current-Changeunit.Currentunit;
            if(Stb.n_current>cycle)
                error('Stb.n_current>cycle');
            elseif(Stb.n_current<1)
                error('Stb.n_current<1');
            end
            end
        end
    end    
    end
    if(flag_CCT==1)
        Exit_thetac=fault.traj.thetac(Stb.n_stbmax,:);
        Exit_omegac=fault.traj.omegac(Stb.n_stbmax,:);
        Exit_theta=fault.traj.theta(Stb.n_stbmax,:);
        Exit_omega=fault.traj.omega(Stb.n_stbmax,:);
    else
        error('No CCT found!');
    end
end
