% =========================================================================
% 教学操作说明：最小梯度点（MGP）搜索或单条轨迹迭代函数
% 使用方法：
%  按以下函数签名调用：function [thetac_MGP,num_Traj,flag_MGP,Normtt,norm_min]=Fun_Cal_MGP(thetac_escape,postfault,preset)
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
function [thetac_MGP,num_Traj,flag_MGP,Normtt,norm_min]=Fun_Cal_MGP(thetac_escape,postfault,preset)
%% Settings
    Tunit=1e-4; % time unit for iteration
    n_itermax=10;    % maximum steps in one iteration procedure
    norm_Tol=1e-5;  % Tolerance set for MGP identification
    n_MGPtraj=0;    % counter for trajectory numbers
    n_MGPtrajmax=1000;  % maximum trajectory numbers
    Yred=postfault.Yred;
    norm_min = 0;
%     ngen=size(Yred,1);
    Normtt=zeros(n_itermax*n_MGPtrajmax,1);  % for observation: all Norm from each trajectory
    f1=evalin('base','f1');
    f2=evalin('base','f2');
%% Initialization
    flag_MGP=0;
    thetac_start=thetac_escape';
    plot(postfault.SEP_delta(2),postfault.SEP_delta(3),'ob','LineWidth',1.5,'MarkerSize',8); hold on;
%% MGP calculation
while(flag_MGP==0)
    [thetac_iter,Normp,no_MGP,flag_MGP]=Fun_Cal_MGP_singletraj(thetac_start',Tunit,n_itermax,norm_Tol,Yred,preset);
    n_MGPtraj=n_MGPtraj+1;  % counter for iteration ++
    for i=1:size(Normp,1)
        noNorm=(n_MGPtraj-1)*n_itermax+i;
        Normtt(noNorm)=Normp(i);
    end
    if(flag_MGP==0)
        figure(f1);
        plot(thetac_iter(:,2),thetac_iter(:,3),'k-','LineWidth',1.5);
        figure(f2);
        plot((n_MGPtraj-1)*n_itermax+1:n_MGPtraj*n_itermax,Normp,'k-','LineWidth',1.5);
        if(n_MGPtraj>1)
            plot((n_MGPtraj-1)*n_itermax:(n_MGPtraj-1)*n_itermax+1,[norm_min Normp(1)],':k','LineWidth',1.5);
        end
        norm_min = Normp(end);
    else
        figure(f1);
        plot(thetac_iter(1:no_MGP,2),thetac_iter(1:no_MGP,3));
        thetac_MGP=thetac_iter(no_MGP,:);
        num_Traj=n_MGPtraj;
        norm_min = Normp(no_MGP);
        figure(f2);
        plot((n_MGPtraj-1)*n_itermax+1:(n_MGPtraj-1)*n_itermax+no_MGP+1,Normp(1:no_MGP+1),'k-','LineWidth',1.5);
        break;
    end
    hold on;
    %% No MGP found in last iteartion process, update start point
    thetac_last=thetac_iter(n_itermax,:)';
    [thetac_update,flag_update]=Fun_Cal_UpdateStartPoint(thetac_last,preset,postfault);
    if(flag_update==1)
        thetac_starthis=thetac_start;
        thetac_start=thetac_update;
%         norm(thetac_start-thetac_starthis)
%         norm(thetac_starthis-postfault.SEP_delta)
        if(norm(thetac_start-thetac_starthis)<1e-3)
            flag_MGP=1;
            thetac_MGP=thetac_update';
            num_Traj=n_MGPtraj;
            fprintf('MGP found since the iteration process reached a repeated status!\n');
            break;
        end
        if(norm(thetac_start-thetac_starthis)>0.5*norm(thetac_starthis-postfault.SEP_delta))
            thetac_start=thetac_last;
            fprintf('No update makes this time\n');
        end
        figure(f1);
        plot([thetac_last(2),thetac_update(2)],[thetac_last(3),thetac_update(3)],':k','LineWidth',1.5);
    else
        thetac_start=thetac_last;
    end
    if(n_MGPtraj>n_MGPtrajmax)
        error('No MGP found in %d times',n_MGPtraj);
    end
end
%% plot the norm of state
%     f2=figure(2);
%     n_count=find(Normtt, 1, 'last' );
%     plot(1:n_count,Normtt(1:n_count));
%     hold on;
%     grid on;


    %%
    clear Normp noNorm i 
