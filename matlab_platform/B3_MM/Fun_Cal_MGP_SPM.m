% =========================================================================
% 教学操作说明：最小梯度点（MGP）搜索或单条轨迹迭代函数
% 使用方法：
%  按以下函数签名调用：function [deltac_MGP,theta_MGP,voltage_MGP,num_Traj,flag_MGP,Normtt,norm_min]=Fun_Cal_MGP_SPM(deltac_escape,theta_escape,voltage_escape,postfault,preset)
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
%% start from escape point, follows stable boundary and find MGP where potential field norm is rather small
%% Input: escape point, postfault network, preset parameters
%% Output: MGP point, numbers of trajectories with each one contains five steps (defined by n_iteamax), flag_MGP (1--success)
function [deltac_MGP,theta_MGP,voltage_MGP,num_Traj,flag_MGP,Normtt,norm_min]=Fun_Cal_MGP_SPM(deltac_escape,theta_escape,voltage_escape,postfault,preset)
%% Settings
    Tunit=1e-3; % time unit for iteration
    n_itermax=10;    % maximum steps in one iteration procedure
    norm_Tol=1e-5;  % Tolerance set for MGP identification
    n_MGPtraj=0;    % counter for trajectory numbers
    n_MGPtrajmax=1000;  % maximum trajectory numbers
    Yfull=postfault.Yfull_mod;
    norm_min = 0;
    Normtt=zeros(n_itermax*n_MGPtrajmax,1);  % for observation: all Norm from each trajectory
     f1=evalin('base','f1');
     f2=evalin('base','f2');
     Basevalue = evalin('base','Basevalue');
%% Initialization
    flag_MGP=0;
    deltac_start=deltac_escape';
    theta_start=theta_escape';
    voltage_start=voltage_escape';
    plot(postfault.SEP_delta(2),postfault.SEP_delta(3),'ob','LineWidth',1.5,'MarkerSize',8); hold on;
%% MGP calculation
while(flag_MGP==0)
    [theta_start,voltage_start,flag_iter,n_iter,err] = Fun_AEiteration_SPM(theta_start,voltage_start,deltac_start,preset,Basevalue,"postfault",1e4,1e-10);
    [deltac_iter,theta_iter,voltage_iter,Normp,no_MGP,flag_MGP]=Fun_Cal_MGP_singletraj_SPM(deltac_start,theta_start,voltage_start,Tunit,n_itermax,norm_Tol,Yfull,preset);
    n_MGPtraj=n_MGPtraj+1;  % counter for iteration ++
    for i=1:size(Normp,1)
        noNorm=(n_MGPtraj-1)*n_itermax+i;
        Normtt(noNorm)=Normp(i);
    end
    if(flag_MGP==0)
        figure(f1);
        plot(deltac_iter(:,2),deltac_iter(:,3),'k-','LineWidth',1.5);
        figure(f2);
        plot((n_MGPtraj-1)*n_itermax+1:n_MGPtraj*n_itermax,Normp,'k-','LineWidth',1.5);
        if(n_MGPtraj>1)
            plot((n_MGPtraj-1)*n_itermax:(n_MGPtraj-1)*n_itermax+1,[norm_min Normp(1)],':k','LineWidth',1.5);
        end
        norm_min = Normp(end);
    else
        figure(f1);
        plot(deltac_iter(1:no_MGP,2),deltac_iter(1:no_MGP,3));
        deltac_MGP=deltac_iter(no_MGP,:);
        theta_MGP=theta_iter(no_MGP,:);
        voltage_MGP= theta_iter(no_MGP,:);
        num_Traj=n_MGPtraj;
        norm_min = Normp(no_MGP);
        figure(f2);
        plot((n_MGPtraj-1)*n_itermax+1:(n_MGPtraj-1)*n_itermax+no_MGP+1,Normp(1:no_MGP+1),'k-','LineWidth',1.5);
        break;
    end
    hold on;
    %% No MGP found in last iteartion process, update start point
    deltac_last=deltac_iter(n_itermax,:)';
    theta_lastpoint=theta_iter(n_itermax,:)';
    voltage_lastpoint=voltage_iter(n_itermax,:)';
    [deltac_update,theta_update,voltage_update,flag_update]=Fun_Cal_UpdateStartPoint_SPM(deltac_last,theta_lastpoint,voltage_lastpoint,preset,postfault);
    if(flag_update==1)
        deltac_starthis=deltac_start;
        deltac_start=deltac_update;
        theta_start=theta_update;
        voltage_start=voltage_update;
        if(norm(deltac_start-deltac_starthis)<1e-3)
            flag_MGP=1;
            deltac_MGP=deltac_update';
            theta_MGP=theta_update';
            voltage_MGP= voltage_update';
            num_Traj=n_MGPtraj;
            fprintf('MGP found since the iteration process reached a repeated status!\n');
            break;
        end
        if(norm(deltac_start-deltac_starthis)>0.5*norm(deltac_starthis-postfault.SEP_delta))
            deltac_start=deltac_last;
            theta_start=theta_lastpoint;
            voltage_start=voltage_lastpoint;
            fprintf('No update makes this time\n');
        end
        figure(f1);
        plot([deltac_last(2),deltac_update(2)],[deltac_last(3),deltac_update(3)],':k','LineWidth',1.5);
    else
        deltac_start=deltac_last;
        theta_start=theta_lastpoint;
        voltage_start=voltage_lastpoint;
    end
    if(n_MGPtraj>n_MGPtrajmax)
        error('No MGP found in %d times',n_MGPtraj);
    end
end
    clear Normp noNorm i 
