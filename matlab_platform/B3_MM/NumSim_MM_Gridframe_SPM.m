% =========================================================================
% 教学操作说明：结构保持模型（SPM）的 DAE 数值轨迹、网络状态和能量绘图脚本。
% 使用方法：
%   新手只通过 run_bcu_beginner.m 的 "spm_numerical" 模式运行；它会先执行
%   Cal_MM_CCT_SPM.m。手动运行前，必须在同一 base workspace 获得 SPM 初始化、
%   CCT 和 postfault.CUEP_delta 结果。
% 参数：
%   无函数形参。Iter.Tfault、Iter.Trecover、Iter.Ttotal、Iter.Tunit 位于下方，
%   单位均为 s，须满足 0 <= Tfault < Trecover <= Ttotal；preset.PathEnergyCal
%   是原始能量计算策略标志，不是物理故障持续时间。
% 返回 / 工作区结果：
%   在 base workspace 生成 DAE 时间域状态、节点相角、节点电压、发电机状态、
%   多类能量项和原始图窗。变量的确切字段以脚本执行后的工作区为准。
% 步骤：
%   1. 读取 SPM 初始化/CCT 结果。2. 计算故障前、故障中和故障后 DAE 轨迹。
%   3. 计算网络电压、节点角度与能量项。4. 绘制机器、网络和能量对比图。
% 单位：
%   Iter 时间变量为 s；机器/网络相角为 rad，角速度为 rad/s；电压、功率和导纳通常为 pu。
% 前置条件：
%   必须已完成同模型 CCT 且存在 CUEP；DAE 轨迹的初始网络状态、节点数和电压
%   维度必须与 Cal_MM_Static_SPM.m 的 case 一致。不要跳过 CCT 直接运行本文件。
% 研究与验证边界：
%   本次仅更新注释，不改变 DAE、积分步长、能量分解或绘图逻辑；计算是否收敛、
%   节点电压是否合理以及能量解释是否成立，仍须以 MATLAB 实际运行和实验记录确认。
% =========================================================================
%close all
%clear
%Cal_MM_Static;
%clear pfdata EMF Yload Case netdata
preset.PathEnergyCal=20; % 0: Ray approximation  n: N-segment trap approximation (-1)--neglect this part
omegab=Basevalue.omegab;
%% 1) 时间设置：SPM 的 fault 区间还要与网络节点保留/删除方式一致。
    Iter.Tfault=20;
    Iter.Trecover=20.2;%20.27;
    Iter.Ttotal=60;
    Iter.Tunit=1e-4;
    T_before = 10;
    T_after = 40;
%% 2) 固定步长轨迹：先对比约简状态，再准备结构保持 DAE 初值。
% prefault
    delta0=prefault.SEP_delta;%zeros(ngen,1)D1D
    omega0=prefault.SEP_omegapu*omegab*ones(ngen,1);
    [theta_pre,omega_pre,thetac_pre,omegacoi_pre,Pe_pre,cycle_pre]=Fun_TrajIter_SRF(Iter.Tfault,Iter.Tunit,prefault.Yred,preset,delta0,omega0,omegab);
% faults
    delta0=theta_pre(cycle_pre,:);
    omega0=omega_pre(cycle_pre,:);
    [theta_fault,omega_fault,thetac_fault,omegacoi_fault,Pe_fault,cycle_fault]=Fun_TrajIter_SRF(Iter.Trecover-Iter.Tfault,Iter.Tunit,fault.Yred,preset,delta0,omega0,omegab);
% postfault
    delta0=theta_fault(cycle_fault,:);
    omega0=omega_fault(cycle_fault,:);
    [theta_post,omega_post,thetac_post,omegacoi_post,Pe_post,cycle_post]=Fun_TrajIter_SRF(Iter.Ttotal-Iter.Trecover,Iter.Tunit,postfault.Yred,preset,delta0,omega0,omegab);
    clear delta0 omega0

%% 3) 约简 ODE 对照：用于与结构保持 DAE 结果比较。
% prefault
    system="prefault";
    delta0=prefault.SEP_delta;
    omega0=prefault.SEP_omegapu*omegab*ones(ngen,1);
    [t_prefault, x_prefault_all] = ode78(@f_timedomain,[0,Iter.Tfault],[delta0; omega0],odeset('RelTol',1e-5));
    deltac_timedomain_prefault=x_prefault_all(:,1:3);
    % faults
    system="fault";
    delta0=x_prefault_all(end,1:3)';
    omega0=x_prefault_all(end,4:6)';
    [t_fault, x_fault_all] = ode78(@f_timedomain,Iter.Tfault:Iter.Tunit:Iter.Trecover,[delta0; omega0],odeset('RelTol',1e-5));
    deltac_timedomain_fault=x_fault_all(:,1:3);

    % postfault
    system="postfault";
    delta0=x_fault_all(end,1:3)';
    omega0=x_fault_all(end,4:6)';
    [t_postfault, x_postfault_all] = ode78(@f_timedomain,[Iter.Trecover,Iter.Ttotal],[delta0; omega0],odeset('RelTol',1e-5));
    clear delta0 omega0
    deltac_timedomain_postfault=x_postfault_all(:,1:3);
    % data collection
    t_timedomain = [t_prefault;t_fault;t_postfault];
    deltac_timedomain = [x_prefault_all(:,1:3); x_fault_all(:,1:3); x_postfault_all(:,1:3)];
    omega_timedomain = [x_prefault_all(:,4:6); x_fault_all(:,4:6); x_postfault_all(:,4:6)]./omegab;
    omegacoi_timedomain = omega_timedomain*preset.m./sum(preset.m);
    omegac_timedomain = omega_timedomain - omegacoi_timedomain*ones(1,3);
    %clear x_prefault_all x_fault_all x_postfault_all

    % bus voltage & angle
    E_pre = preset.Epu'.*(cos(deltac_timedomain_prefault)+1i*sin(deltac_timedomain_prefault));
    Vnet_pre = -(inv(prefault.Yrr)*prefault.Yrn*E_pre.').';
    Vmnet_pre = abs(Vnet_pre);
    Vpnet_pre = angle(Vnet_pre);
    
    E_fault = preset.Epu'.*(cos(deltac_timedomain_fault)+1i*sin(deltac_timedomain_fault));
    Vnet_fault = -(inv(fault.Yrr)*fault.Yrn*E_fault.').';
    temp_Vmnet_fault = abs(Vnet_fault);
    Vmnet_fault=zeros(size(temp_Vmnet_fault,1),6);
    temp_Vpnet_fault = angle(Vnet_fault);
    Vpnet_fault= zeros(size(temp_Vpnet_fault,1),6);
    for i=1:6
       if (i+3)<fault.faultbus
           Vmnet_fault(:,i)=temp_Vmnet_fault(:,i);
           Vpnet_fault(:,i)=temp_Vpnet_fault(:,i);
       elseif (i+3)>fault.faultbus
           Vmnet_fault(:,i)=temp_Vmnet_fault(:,i-1);
           Vpnet_fault(:,i)=temp_Vpnet_fault(:,i-1);
       end
    end
    
    E_post = preset.Epu'.*(cos(deltac_timedomain_postfault)+1i*sin(deltac_timedomain_postfault));
    Vnet_post = -(inv(postfault.Yrr)*postfault.Yrn*E_post.').';
    Vmnet_post = abs(Vnet_post);
    Vpnet_post = angle(Vnet_post);
    
    Vmnet = [Vmnet_pre;Vmnet_fault;Vmnet_post];
    Vpnet = [Vpnet_pre;Vpnet_fault;Vpnet_post];
%% 4) 结构保持 DAE：Mass 矩阵把网络代数约束与发电机微分方程组合起来。
M = diag([ones(3,1); ones(12,1)*1e-15; ones(3,1)]);
options = odeset('Mass',M,'RelTol',1e-10,'AbsTol',[1e-8*ones(1,3),1e-12*ones(1,12),1e-8*ones(1,3)]);
% prefault
    system="prefault";
    delta0=prefault.SEP_delta;
    omega0=prefault.SEP_omegapu*omegab*ones(ngen,1);
    delta_net0=prefault.net_delta;
    voltage_net0=prefault.net_voltage;
    x0=[delta0; delta_net0; voltage_net0; omega0];
    [t_prefault, x_prefault_all] = ode15s(@(t,x)f_timedomain_DAE(t, x, system),0:Iter.Tunit:Iter.Tfault,x0,options);
%fault
    system="fault1";
    delta0=x_prefault_all(end,1:3)';
    omega0=x_prefault_all(end,16:18)';
    delta_net0=x_prefault_all(end,4:9)';
    voltage_net0=x_prefault_all(end,10:15)';
    delta_net0(fault.faultbus-ngen)=[];
    voltage_net0(fault.faultbus-ngen)=[];
    delta_net0(6)=0;
    voltage_net0(6)=0;
    [delta_net_s,V_net_s,flag_iter,n_iter,err] = Fun_AEiteration_SPM(delta_net0,voltage_net0,delta0,preset,Basevalue,system,1e5,1e-10);
    [t_fault, x_fault_all] = ode15s(@(t, x) f_timedomain_DAE(t, x, system),Iter.Tfault:Iter.Tunit:Iter.Trecover,[delta0; delta_net_s; V_net_s; omega0],options);
    fault_delta_net=x_fault_all(2:end,4:8);
    fault_delta_net=[fault_delta_net(:,1:(fault.faultbus-ngen-1)) zeros(size(fault_delta_net,1),1) fault_delta_net(:,(fault.faultbus-ngen):end)];
    fault_voltage_net=x_fault_all(2:end,10:14);
    fault_voltage_net=[fault_voltage_net(:,1:(fault.faultbus-ngen-1)) zeros(size(fault_voltage_net,1),1) fault_voltage_net(:,(fault.faultbus-ngen):end)];
    fault_delta_gen = x_fault_all(2:end,1:3);
%% fault-clear
    system = "postfault";
    delta_net_faultclear = zeros(size(fault_delta_gen,1),nbus-ngen);
    voltage_net_faultclear = zeros(size(fault_delta_gen,1),nbus-ngen);
    delta_net_faultclear2 = zeros(size(fault_delta_gen,1),nbus-ngen);
    voltage_net_faultclear2 = zeros(size(fault_delta_gen,1),nbus-ngen);
    temp_ini = [zeros(6,1);ones(6,1)];
    for i = 1:(size(fault_delta_gen,1))
       %[delta_net_temp,voltage_temp,flag_iter,n_iter,err] = Fun_AEiteration_SPM(zeros(6,1),ones(6,1),fault_delta_gen(i,:)',preset,Basevalue,system,1e4,1e-10);
       %delta_net_faultclear(i,:) = delta_net_temp';
       %voltage_net_faultclear(i,:) = voltage_temp';
       [Results_fsolve,fval,exitflag,output]=fsolve(@(x)Fun_AEfslove_SPM(x,fault_delta_gen(i,:)',preset,system),temp_ini,optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off','TolX',1e-5));
       if exitflag<=0
           error('cannot find fault-clear state! \n');
       end
       delta_net_faultclear(i,:) = Results_fsolve(1:6)';
       voltage_net_faultclear(i,:) = Results_fsolve(7:12)';
       temp_ini = Results_fsolve;
    end

%% postfault
    system="postfault";
    delta0=x_fault_all(end,1:3)';
    omega0=x_fault_all(end,16:18)';
    [delta_net_s,Voltage_net_s,flag_iter,n_iter,err] = Fun_AEiteration_SPM(temp_ini(1:6),temp_ini(7:12),delta0,preset,Basevalue,system,1e4,1e-10);
    [Results_fsolve,fval,exitflag,output]=fsolve(@(x)Fun_AEfslove_SPM(x,delta0,preset,system),temp_ini,optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off','TolX',1e-5));
    
    [t_postfault, x_postfault_all] = ode15s(@(t,x)f_timedomain_DAE(t,x,system),Iter.Trecover:Iter.Tunit:Iter.Ttotal,[delta0; delta_net_s; Voltage_net_s; omega0],options);
    %clear delta0 omega0

% data collection
     t_timedomain_DAE = [t_prefault(2:end);t_fault(2:end);t_postfault(2:end)];
     deltac_timedomain_DAE  = [x_prefault_all(2:end,1:3); fault_delta_gen; x_postfault_all(2:end,1:3)];
     omega_timedomain_DAE  = [x_prefault_all(2:end,16:18); x_fault_all(2:end,16:18); x_postfault_all(2:end,16:18)]./omegab;
     omegacoi_timedomain_DAE  = omega_timedomain_DAE *preset.m./sum(preset.m);
     omegac_timedomain_DAE  = omega_timedomain_DAE  - omegacoi_timedomain_DAE *ones(1,3);
     delta_net_DAE = [x_prefault_all(2:end,4:9); fault_delta_net; x_postfault_all(2:end,4:9)];
     voltage_net_DAE = [x_prefault_all(2:end,10:15); fault_voltage_net; x_postfault_all(2:end,10:15)];

     delta_net_clear_DAE = [x_prefault_all(2:end,4:9); delta_net_faultclear; x_postfault_all(2:end,4:9)];
     voltage_net_clear_DAE = [x_prefault_all(2:end,10:15); voltage_net_faultclear; x_postfault_all(2:end,10:15)];

%% dae based Runge-Kutta

% system = "postfault";
% delta0=x_fault_all(end,1:3)';
% omega0=x_fault_all(end,16:18)';
% [delta_net_s,Voltage_net_s,flag_iter,n_iter,err] = Fun_AEiteration_SPM(temp_ini(1:6),temp_ini(7:12),delta0,preset,Basevalue,system,1e4,1e-10);
% [delta_RK4,omega_RK4,deltac_RK4,omegacoi_RK4,theta_net_RK4,voltage_net_RK4,cycle]=Fun_TrajIter_SPM(Iter.Ttotal-Iter.Trecover,Iter.Tunit,system,preset,delta0,omega0,delta_net_s,Voltage_net_s,Basevalue);
% omegac_RK4=(omega_RK4-omegacoi_RK4)./omegab;
% omegacoi_RK4=omegacoi_RK4./omegab;
% 
% Tpre=Iter.Tfault;
% Tfault=Iter.Trecover-Iter.Tfault;
% Tpost=Iter.Ttotal-Iter.Trecover;
% Tunit=Iter.Tunit;
% Ttotal=Iter.Ttotal;
% TM_pre=0:Tunit:Tpre-Tunit;
% TM_fault=Tpre:Tunit:Tpre+Tfault-Tunit;
% TM_post=round((Tpre+Tfault)/Tunit)*Tunit:Tunit:(round((Ttotal)/Tunit)-1)*Tunit;

%% plot for DAE
% deltac
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE,deltac_timedomain_DAE(:,1),'linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
    plot(t_timedomain_DAE,deltac_timedomain_DAE(:,2),'linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
    plot(t_timedomain_DAE,deltac_timedomain_DAE(:,3),'linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
% 
%         plot(IterData.Tout(n_start:n_end),IterData.thetac((n_start:n_end),1),'LineStyle','--','linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
%         plot(IterData.Tout(n_start:n_end),IterData.thetac((n_start:n_end),2),'LineStyle','--','linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
%         plot(IterData.Tout(n_start:n_end),IterData.thetac((n_start:n_end),3),'LineStyle','--','linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
%        

%      plot(TM_post,delta_RK4(:,1),'LineStyle',':','linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
%      plot(TM_post,delta_RK4(:,2),'LineStyle',':','linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
%      plot(TM_post,delta_RK4(:,3),'LineStyle',':','linewidth',2,'color',[0/255 175/255 0/255]);    hold on;

    grid on;    grid minor;
    ylim([-2,3]);
    yl=ylim;
    ymin=yl(1,1);
    ymax=yl(1,2);
    xlim([20-T_before 20+T_after]);
    thetarange=[ymin,ymin,ymax,ymax];
    xl=xlim;
    xmin=xl(1,1);
    xmax=xl(1,2);
    ylab_x=xmin-(xmax-xmin)/15;
    ylab_y=(ymax+ymin)/2;
    xlab_x=(xmax+xmin)/2;
    xlab_y=ymin-(ymax-ymin)/15;
    ylim([ymin ymax]);
    trange=[Iter.Tfault,Iter.Trecover,Iter.Trecover,Iter.Tfault];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    axis([xmin xmax ymin ymax ]);
    % axis font
    ax=gca;
    ax.FontName='Arial';
    ax.FontSize=14;
    xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
    ylabel('\delta_c(rad)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);
%% omega
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);

    plot(t_timedomain_DAE,omegac_timedomain_DAE(:,1),'linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
    plot(t_timedomain_DAE,omegac_timedomain_DAE(:,2),'linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
    plot(t_timedomain_DAE,omegac_timedomain_DAE(:,3),'linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
    plot(t_timedomain_DAE,omegacoi_timedomain_DAE-1,'linewidth',2,'color',[0/255 0/255 0/255]);    hold on;


%     plot(TM_post,omegac_RK4(:,1),'LineStyle',':','linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
%     plot(TM_post,omegac_RK4(:,2),'LineStyle',':','linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
%     plot(TM_post,omegac_RK4(:,3),'LineStyle',':','linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
%     plot(TM_post,omegacoi_RK4-1,'LineStyle',':','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;

    xlim([20-T_before 20+T_after]);
    grid on;    grid minor;
    ylim([-0.03,0.05]);

    yl=ylim;
    ymin=yl(1,1);
    ymax=yl(1,2);
    xl=xlim;
    xmin=xl(1,1);
    xmax=xl(1,2);
    ylab_x=xmin-(xmax-xmin)/12;
    ylab_y=(ymax+ymin)/2;
    xlab_x=(xmax+xmin)/2;
    xlab_y=ymin-(ymax-ymin)/15;
    % during-fault area identification
    trange=[Iter.Tfault,Iter.Trecover,Iter.Trecover,Iter.Tfault];   thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
    ylim([ymin,ymax]);
    % axis font
    ax=gca;
    ax.FontName='Arial';
    ax.FontSize=14;
    xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
    ylabel('\omega(pu)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);
    
    %% Bus angle
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);

   % plot(t_timedomain,deltac_timedomain(:,1),'linewidth',2);    hold on;
   % plot(t_timedomain,deltac_timedomain(:,2),'linewidth',2);    hold on;
   % plot(t_timedomain,deltac_timedomain(:,3),'linewidth',2);    hold on;

    plot(t_timedomain_DAE,delta_net_DAE(:,1),'linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
    plot(t_timedomain_DAE,delta_net_DAE(:,2),'linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
    plot(t_timedomain_DAE,delta_net_DAE(:,3),'linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
    plot(t_timedomain_DAE,delta_net_DAE(:,4),'linewidth',2,'color',[0.4660 0.6740 0.1880]);    hold on;
    plot(t_timedomain_DAE,delta_net_DAE(:,5),'linewidth',2,'color',[0.3010 0.7450 0.9330]);    hold on;
    plot(t_timedomain_DAE,delta_net_DAE(:,6),'linewidth',2,'color',[0.6350 0.0780 0.1840]);    hold on;

    plot(t_fault(2:end),delta_net_faultclear(:,1),'LineStyle',':','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
    plot(t_fault(2:end),delta_net_faultclear(:,2),'LineStyle',':','linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
    plot(t_fault(2:end),delta_net_faultclear(:,3),'LineStyle',':','linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
    plot(t_fault(2:end),delta_net_faultclear(:,4),'LineStyle',':','linewidth',2,'color',[0.4660 0.6740 0.1880]);    hold on;
    plot(t_fault(2:end),delta_net_faultclear(:,5),'LineStyle',':','linewidth',2,'color',[0.3010 0.7450 0.9330]);    hold on;
    plot(t_fault(2:end),delta_net_faultclear(:,6),'LineStyle',':','linewidth',2,'color',[0.6350 0.0780 0.1840]);    hold on;


%     plot(TM_post,theta_net_RK4(:,1),'LineStyle',':','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
%     plot(TM_post,theta_net_RK4(:,2),'LineStyle',':','linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
%     plot(TM_post,theta_net_RK4(:,3),'LineStyle',':','linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
%     plot(TM_post,theta_net_RK4(:,4),'LineStyle',':','linewidth',2,'color',[0.4660 0.6740 0.1880]);    hold on;
%     plot(TM_post,theta_net_RK4(:,5),'LineStyle',':','linewidth',2,'color',[0.3010 0.7450 0.9330]);    hold on;
%     plot(TM_post,theta_net_RK4(:,6),'LineStyle',':','linewidth',2,'color',[0.6350 0.0780 0.1840]);    hold on;


    
    xlim([20-T_before 20+T_after]);
    grid on;    grid minor;
    ylim([-2,3]);
    % xlabel ylabel position
    yl=ylim;
    ymin=yl(1,1);
    ymax=yl(1,2);
    xl=xlim;
    xmin=xl(1,1);
    xmax=xl(1,2);
    ylab_x=xmin-(xmax-xmin)/15;
    ylab_y=(ymax+ymin)/2;
    xlab_x=(xmax+xmin)/2;
    xlab_y=ymin-(ymax-ymin)/15;
    % during-fault area identification
    trange=[Iter.Tfault,Iter.Trecover,Iter.Trecover,Iter.Tfault];   thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
    ylim([ymin,ymax]);
    % axis font
    ax=gca;
    ax.FontName='Arial';
    ax.FontSize=14;
    legendEntries = cell(1, 6);
    for i = 1:6
        legendEntries{i} = ['Bus' num2str(i+3)];
    end
    legend(legendEntries{:});
    xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
    ylabel('\theta_c(rad)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);

    %% Bus voltage
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    
  %  plot(t_timedomain,E(1)*ones(size(t_timedomain,1),1),'linewidth',2);    hold on;
  %  plot(t_timedomain,E(2)*ones(size(t_timedomain,1),1),'linewidth',2);    hold on;
  %  plot(t_timedomain,E(3)*ones(size(t_timedomain,1),1),'linewidth',2);    hold on;

    plot(t_timedomain_DAE,voltage_net_DAE(:,1),'linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
    plot(t_timedomain_DAE,voltage_net_DAE(:,2),'linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
    plot(t_timedomain_DAE,voltage_net_DAE(:,3),'linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
    plot(t_timedomain_DAE,voltage_net_DAE(:,4),'linewidth',2,'color',[0.4660 0.6740 0.1880]);    hold on;
    plot(t_timedomain_DAE,voltage_net_DAE(:,5),'linewidth',2,'color',[0.3010 0.7450 0.9330]);    hold on;
    plot(t_timedomain_DAE,voltage_net_DAE(:,6),'linewidth',2,'color',[0.6350 0.0780 0.1840]);    hold on;

    plot(t_fault(2:end),voltage_net_faultclear(:,1),'LineStyle',':','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
    plot(t_fault(2:end),voltage_net_faultclear(:,2),'LineStyle',':','linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
    plot(t_fault(2:end),voltage_net_faultclear(:,3),'LineStyle',':','linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
    plot(t_fault(2:end),voltage_net_faultclear(:,4),'LineStyle',':','linewidth',2,'color',[0.4660 0.6740 0.1880]);    hold on;
    plot(t_fault(2:end),voltage_net_faultclear(:,5),'LineStyle',':','linewidth',2,'color',[0.3010 0.7450 0.9330]);    hold on;
    plot(t_fault(2:end),voltage_net_faultclear(:,6),'LineStyle',':','linewidth',2,'color',[0.6350 0.0780 0.1840]);    hold on;


%     plot(TM_post,voltage_net_RK4(:,1),'LineStyle',':','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
%     plot(TM_post,voltage_net_RK4(:,2),'LineStyle',':','linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
%     plot(TM_post,voltage_net_RK4(:,3),'LineStyle',':','linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
%     plot(TM_post,voltage_net_RK4(:,4),'LineStyle',':','linewidth',2,'color',[0.4660 0.6740 0.1880]);    hold on;
%     plot(TM_post,voltage_net_RK4(:,5),'LineStyle',':','linewidth',2,'color',[0.3010 0.7450 0.9330]);    hold on;
%     plot(TM_post,voltage_net_RK4(:,6),'LineStyle',':','linewidth',2,'color',[0.6350 0.0780 0.1840]);    hold on;
%     
    xlim([20-T_before 20+T_after]);
    grid on;    grid minor;
    ylim([-0.1,1.2]);
    % xlabel ylabel position
    yl=ylim;
    ymin=yl(1,1);
    ymax=yl(1,2);
    xl=xlim;
    xmin=xl(1,1);
    xmax=xl(1,2);
    ylab_x=xmin-(xmax-xmin)/15;
    ylab_y=(ymax+ymin)/2;
    xlab_x=(xmax+xmin)/2;
    xlab_y=ymin-(ymax-ymin)/15;
    % during-fault area identification
    trange=[Iter.Tfault,Iter.Trecover,Iter.Trecover,Iter.Tfault];   thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
    ylim([ymin,ymax]);
    % axis font
    ax=gca;
    ax.FontName='Arial';
    ax.FontSize=14;
    legendEntries = cell(1, 6);
    for i = 1:6
        legendEntries{i} = ['Bus' num2str(i+3)];
    end
    legend(legendEntries{:});
    xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
    ylabel('voltage(p.u.)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);

%% 5) 数据拼接：保留发电机角度、网络角度、电压和 COI 量的对应索引。
    Tpre=Iter.Tfault;
    Tfault=Iter.Trecover-Iter.Tfault;
    Tpost=Iter.Ttotal-Iter.Trecover;
    Tunit=Iter.Tunit;
    Ttotal=Iter.Ttotal;
    TM_pre=0:Tunit:Tpre-Tunit;
    TM_fault=Tpre:Tunit:Tpre+Tfault-Tunit;
    TM_post=round((Tpre+Tfault)/Tunit)*Tunit:Tunit:(round((Ttotal)/Tunit)-1)*Tunit;
    clear Tpre Tfault Tpost Tunit
    n_tt=cycle_pre+cycle_fault+cycle_post;
    IterData.omega=zeros(n_tt,ngen);
    IterData.omegacoi=zeros(n_tt,1);
    IterData.omegac=zeros(n_tt,ngen);
    IterData.theta=zeros(n_tt,ngen);
    IterData.thetac=zeros(n_tt,ngen);
    IterData.Tout=zeros(n_tt,1);
    IterData.Pe=zeros(n_tt,ngen);
    IterData.theta(1:cycle_pre,:)=theta_pre;
    IterData.thetac(1:cycle_pre,:)=thetac_pre;
    IterData.omega(1:cycle_pre,:)=omega_pre/omegab;
    IterData.omegacoi(1:cycle_pre)=omegacoi_pre/omegab;
    for i=1:ngen
        IterData.omegac(1:cycle_pre,i)=omega_pre(:,i)/omegab-omegacoi_pre/omegab;
    end
    IterData.Tout(1:cycle_pre)=TM_pre;
    IterData.Pe(1:cycle_pre,:)=Pe_pre;

    IterData.theta(cycle_pre+1:cycle_pre+cycle_fault,:)=theta_fault;
    IterData.thetac(cycle_pre+1:cycle_pre+cycle_fault,:)=thetac_fault;
    IterData.omega(cycle_pre+1:cycle_pre+cycle_fault,:)=omega_fault/omegab;
    IterData.omegacoi(cycle_pre+1:cycle_pre+cycle_fault)=omegacoi_fault/omegab;
    for i=1:ngen
        IterData.omegac(cycle_pre+1:cycle_pre+cycle_fault,i)=omega_fault(:,i)/omegab-omegacoi_fault/omegab;
    end
    IterData.Tout(cycle_pre+1:cycle_pre+cycle_fault)=TM_fault;
    IterData.Pe(cycle_pre+1:cycle_pre+cycle_fault,:)=Pe_fault;

    IterData.theta(cycle_pre+cycle_fault+1:n_tt,:)=theta_post;
    IterData.thetac(cycle_pre+cycle_fault+1:n_tt,:)=thetac_post;
    IterData.omega(cycle_pre+cycle_fault+1:n_tt,:)=omega_post/omegab;
    IterData.omegacoi(cycle_pre+cycle_fault+1:n_tt)=omegacoi_post/omegab;
    for i=1:ngen
        IterData.omegac(cycle_pre+cycle_fault+1:n_tt,i)=omega_post(:,i)/omegab-omegacoi_post/omegab;
    end
    IterData.omegacd=IterData.omegac*preset.d/sum(preset.d);
    IterData.deltacd=IterData.thetac*preset.d/sum(preset.d);
    
    IterData.Tout(cycle_pre+cycle_fault+1:n_tt)=TM_post;
    IterData.Pe(cycle_pre+cycle_fault+1:n_tt,:)=Pe_post;
    clear n_tt


    %% Plot original
    % thetac
%         figure;
%         set(gca,'position',[0.115,0.12,0.815,0.84]);
%         set(gcf,'position',[60 200 600 450]);

        n_start=cycle_pre-fix(T_before/Iter.Tunit)+1;
        n_end=cycle_pre+fix(T_after/Iter.Tunit);
        
%         plot(IterData.Tout(n_start:n_end),IterData.thetac((n_start:n_end),1),'linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
%         plot(IterData.Tout(n_start:n_end),IterData.thetac((n_start:n_end),2),'linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
%         plot(IterData.Tout(n_start:n_end),IterData.thetac((n_start:n_end),3),'linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
%         plot(IterData.Tout(n_start:n_end),IterData.deltacd((n_start:n_end),1),'linewidth',2,'color',[200/255 200/255 200/255]);    hold on;
% 
%         %ode time domain
% %         plot(t_timedomain,deltac_timedomain(:,1),'LineStyle',':','linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
% %         plot(t_timedomain,deltac_timedomain(:,2),'LineStyle',':','linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
% %         plot(t_timedomain,deltac_timedomain(:,3),'LineStyle',':','linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
%      plot(TM_post,delta_RK4(:,1),'LineStyle',':','linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
%      plot(TM_post,delta_RK4(:,2),'LineStyle',':','linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
%      plot(TM_post,delta_RK4(:,3),'LineStyle',':','linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
% 
%         xlim([20-T_before 20+T_after]);
%         grid on;    grid minor;
%         ylim([-2,4]);
%         % xlabel ylabel position
%         yl=ylim;
%         ymin=yl(1,1);
%         ymax=yl(1,2);
%         xl=xlim;
%         xmin=xl(1,1);
%         xmax=xl(1,2);
%         ylab_x=xmin-(xmax-xmin)/15;
%         ylab_y=(ymax+ymin)/2;
%         xlab_x=(xmax+xmin)/2;
%         xlab_y=ymin-(ymax-ymin)/15;
%         % during-fault area identification
%         trange=[Iter.Tfault,Iter.Trecover,Iter.Trecover,Iter.Tfault];   thetarange=[ymin,ymin,ymax,ymax];
%         fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
%         ylim([ymin,ymax]);
%         % axis font
%         ax=gca;
%         ax.FontName='Arial';
%         ax.FontSize=14;
%         xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
%         ylabel('\delta_c(rad)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);


%         figure;
%         set(gca,'position',[0.115,0.12,0.815,0.84]);
%         set(gcf,'position',[60 200 600 450]);
%         plot(IterData.Tout(n_start:n_end),IterData.omegac((n_start:n_end),1),'linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
%         plot(IterData.Tout(n_start:n_end),IterData.omegac((n_start:n_end),2),'linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
%         plot(IterData.Tout(n_start:n_end),IterData.omegac((n_start:n_end),3),'linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
%         
%         plot(IterData.Tout(n_start:n_end),IterData.omegacoi((n_start:n_end),1)-1,'linewidth',2,'color',[0/255 0/255 0/255]);    hold on;
%         plot(IterData.Tout(n_start:n_end),IterData.omegacd((n_start:n_end),1),'linewidth',2,'color',[200/255 200/255 200/255]);    hold on;
% 
% %         plot(t_timedomain,omegac_timedomain(:,1),'LineStyle',':','linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
% %         plot(t_timedomain,omegac_timedomain(:,2),'LineStyle',':','linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
% %         plot(t_timedomain,omegac_timedomain(:,3),'LineStyle',':','linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
% %         plot(t_timedomain,omegacoi_timedomain-1,'LineStyle',':','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;
%         plot(TM_post,omegac_RK4(:,1),'LineStyle',':','linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
%         plot(TM_post,omegac_RK4(:,2),'LineStyle',':','linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
%         plot(TM_post,omegac_RK4(:,3),'LineStyle',':','linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
%         plot(TM_post,omegacoi_RK4-1,'LineStyle',':','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;
% 
%         xlim([20-T_before 20+T_after]);
%         grid on;    grid minor;
%         ylim([-0.03,0.05]);
% 
%         yl=ylim;
%         ymin=yl(1,1);
%         ymax=yl(1,2);
%         xl=xlim;
%         xmin=xl(1,1);
%         xmax=xl(1,2);
%         ylab_x=xmin-(xmax-xmin)/12;
%         ylab_y=(ymax+ymin)/2;
%         xlab_x=(xmax+xmin)/2;
%         xlab_y=ymin-(ymax-ymin)/15;
%         % during-fault area identification
%         trange=[Iter.Tfault,Iter.Trecover,Iter.Trecover,Iter.Tfault];   thetarange=[ymin,ymin,ymax,ymax];
%         fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
%         ylim([ymin,ymax]);
%         % axis font
%         ax=gca;
%         ax.FontName='Arial';
%         ax.FontSize=14;
%         xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
%         ylabel('\omega(pu)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);
%     
%         %Bus angle
%         figure;
%         set(gca,'position',[0.115,0.12,0.815,0.84]);
%         set(gcf,'position',[60 200 600 450]);
% 
%        % plot(t_timedomain,deltac_timedomain(:,1),'linewidth',2);    hold on;
%        % plot(t_timedomain,deltac_timedomain(:,2),'linewidth',2);    hold on;
%        % plot(t_timedomain,deltac_timedomain(:,3),'linewidth',2);    hold on;
% 
%         plot(t_timedomain,Vpnet(:,1),'linewidth',2);    hold on;
%         plot(t_timedomain,Vpnet(:,2),'linewidth',2);    hold on;
%         plot(t_timedomain,Vpnet(:,3),'linewidth',2);    hold on;
%         plot(t_timedomain,Vpnet(:,4),'linewidth',2);    hold on;
%         plot(t_timedomain,Vpnet(:,5),'linewidth',2);    hold on;
%         plot(t_timedomain,Vpnet(:,6),'linewidth',2);    hold on;
% 
%         plot(TM_post,theta_net_RK4(:,1),'LineStyle',':','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
%         plot(TM_post,theta_net_RK4(:,2),'LineStyle',':','linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
%         plot(TM_post,theta_net_RK4(:,3),'LineStyle',':','linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
%         plot(TM_post,theta_net_RK4(:,4),'LineStyle',':','linewidth',2,'color',[0.4660 0.6740 0.1880]);    hold on;
%         plot(TM_post,theta_net_RK4(:,5),'LineStyle',':','linewidth',2,'color',[0.3010 0.7450 0.9330]);    hold on;
%         plot(TM_post,theta_net_RK4(:,6),'LineStyle',':','linewidth',2,'color',[0.6350 0.0780 0.1840]);    hold on;
% 
%         
%         xlim([20-T_before 20+T_after]);
%         grid on;    grid minor;
%         ylim([-2,4]);
%         % xlabel ylabel position
%         yl=ylim;
%         ymin=yl(1,1);
%         ymax=yl(1,2);
%         xl=xlim;
%         xmin=xl(1,1);
%         xmax=xl(1,2);
%         ylab_x=xmin-(xmax-xmin)/15;
%         ylab_y=(ymax+ymin)/2;
%         xlab_x=(xmax+xmin)/2;
%         xlab_y=ymin-(ymax-ymin)/15;
%         % during-fault area identification
%         trange=[Iter.Tfault,Iter.Trecover,Iter.Trecover,Iter.Tfault];   thetarange=[ymin,ymin,ymax,ymax];
%         fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
%         ylim([ymin,ymax]);
%         % axis font
%         ax=gca;
%         ax.FontName='Arial';
%         ax.FontSize=14;
%         legendEntries = cell(1, 6);
%         for i = 1:6
%             legendEntries{i} = ['Bus' num2str(i+3)];
%         end
%         legend(legendEntries{:});
%         xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
%         ylabel('\delta_c(rad)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);
%     
%         %Bus voltage
%         figure;
%         set(gca,'position',[0.115,0.12,0.815,0.84]);
%         set(gcf,'position',[60 200 600 450]);
%         
%       %  plot(t_timedomain,E(1)*ones(size(t_timedomain,1),1),'linewidth',2);    hold on;
%       %  plot(t_timedomain,E(2)*ones(size(t_timedomain,1),1),'linewidth',2);    hold on;
%       %  plot(t_timedomain,E(3)*ones(size(t_timedomain,1),1),'linewidth',2);    hold on;
% 
%         plot(t_timedomain,Vmnet(:,1),'linewidth',2);    hold on;
%         plot(t_timedomain,Vmnet(:,2),'linewidth',2);    hold on;
%         plot(t_timedomain,Vmnet(:,3),'linewidth',2);    hold on;
%         plot(t_timedomain,Vmnet(:,4),'linewidth',2);    hold on;
%         plot(t_timedomain,Vmnet(:,5),'linewidth',2);    hold on;
%         plot(t_timedomain,Vmnet(:,6),'linewidth',2);    hold on;
%                 
%         plot(TM_post,voltage_net_RK4(:,1),'LineStyle',':','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
%         plot(TM_post,voltage_net_RK4(:,2),'LineStyle',':','linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
%         plot(TM_post,voltage_net_RK4(:,3),'LineStyle',':','linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
%         plot(TM_post,voltage_net_RK4(:,4),'LineStyle',':','linewidth',2,'color',[0.4660 0.6740 0.1880]);    hold on;
%         plot(TM_post,voltage_net_RK4(:,5),'LineStyle',':','linewidth',2,'color',[0.3010 0.7450 0.9330]);    hold on;
%         plot(TM_post,voltage_net_RK4(:,6),'LineStyle',':','linewidth',2,'color',[0.6350 0.0780 0.1840]);    hold on;
%         
%         xlim([20-T_before 20+T_after]);
%         grid on;    grid minor;
%         ylim([-0.1,1.2]);
%         % xlabel ylabel position
%         yl=ylim;
%         ymin=yl(1,1);
%         ymax=yl(1,2);
%         xl=xlim;
%         xmin=xl(1,1);
%         xmax=xl(1,2);
%         ylab_x=xmin-(xmax-xmin)/15;
%         ylab_y=(ymax+ymin)/2;
%         xlab_x=(xmax+xmin)/2;
%         xlab_y=ymin-(ymax-ymin)/15;
%         % during-fault area identification
%         trange=[Iter.Tfault,Iter.Trecover,Iter.Trecover,Iter.Tfault];   thetarange=[ymin,ymin,ymax,ymax];
%         fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
%         ylim([ymin,ymax]);
%         % axis font
%         ax=gca;
%         ax.FontName='Arial';
%         ax.FontSize=14;
%         legendEntries = cell(1, 6);
%         for i = 1:6
%             legendEntries{i} = ['Bus' num2str(i+3)];
%         end
%         legend(legendEntries{:});
%         xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
%         ylabel('\delta_c(rad)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);


%% Calculate Energy
     no_duration = n_end - n_start + 1;
     Pm=preset.Pmpu;
     E=preset.Epu;
     m=preset.m;
     d=preset.d;
     Yfull_post=postfault.Yfull_mod;
     G_post=real(Yfull_post);
     B_post=imag(Yfull_post);
    % Potential energy
        [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3),Ep_tmp(4),Ep_tmp(5)]=Fun_Cal_PotentialEnergy_SPM(preset,postfault,deltac_timedomain_DAE(n_start,:)',delta_net_clear_DAE(n_start,:)',voltage_net_clear_DAE(n_start,:)');
        Ep_start=sum(Ep_tmp);
        Ep_loadloss_start=Ep_tmp(3);
        Ep_mag_start=Ep_tmp(2);
        Ep_mag_load_start = Ep_tmp(5);
        Ep_lossy_start = Ep_tmp(4);
        clear Ep_tmp
        Ep=zeros(no_duration,1);
        Ep_loadloss=zeros(no_duration,1);
        Ep_lossy=zeros(no_duration,1);
        Ep_mag=zeros(no_duration,1);
        Ep_mag_load=zeros(no_duration,1);
        for tm=1:size(Ep,1)
            [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3),Ep_tmp(4),Ep_tmp(5)]=Fun_Cal_PotentialEnergy_SPM(preset,postfault,deltac_timedomain_DAE(n_start-1+tm,:)',delta_net_clear_DAE(n_start-1+tm,:)',voltage_net_clear_DAE(n_start-1+tm,:)');
            Ep_loadloss(tm)=Ep_tmp(3);
            Ep_mag(tm)=Ep_tmp(2);
            Ep_lossy(tm)=Ep_tmp(4);
            Ep_mag_load(tm)=Ep_tmp(5);
            Ep(tm)=sum(Ep_tmp);   clear Ep_tmp
        end

     %% Kinetic energy
        Ek0=zeros(ngen,1);
        Eke=zeros(ngen,1);
        for i=1:ngen
            Ek0(i)=0.5*m(i)*(omegac_timedomain_DAE(n_start,i)*omegab)^2;
        end
        Ek_start=sum(Ek0);
    % Damping energy in iteration
        % for Potential energy error display
        Ed_non_record=zeros(no_duration,1);
        Ed_un_record=zeros(no_duration,1);
        Ep_lossy_iter=zeros(no_duration,1);
        Ep_loadloss_iter=zeros(no_duration,1);
        Ep_magload_iter=zeros(no_duration,1);
        Ep_mag_iter=zeros(no_duration,1);
        Ep_iter_record=zeros(no_duration,1);
        Ep_iter_record_withoutnet=zeros(no_duration,1);
        P_un_record=zeros(no_duration,1);
        P_non_record=zeros(no_duration,1);
        Ep_iter_record(1)=Ep_start;
        Ep_iter_record_withoutnet(1)=Ep_start;
        Ep_lossy_iter(1)=Ep_lossy_start;
        Ep_loadloss_iter(1) = Ep_loadloss_start;
        Ep_magload_iter(1) = Ep_mag_load_start;
        Ep_mag_iter(1)=Ep_mag_start;
        Ek=zeros(no_duration,1);
        Ek(1)=Ek_start;
        Err_step=zeros(no_duration,1);
        Err_step_exp=zeros(no_duration,1);

        Pe_pre=zeros(ngen,1);
        Pnet_pre=zeros(nbus-ngen,1);
        Qnet_pre=zeros(nbus-ngen,1);
        Pe_lossy_pre=zeros(ngen,1);  % lossy dispattive energy
        Pnet_lossy_pre=zeros(nbus-ngen,1);  % lossy dispattive energy
        Qnet_lossy_pre=zeros(nbus-ngen,1);  % lossy dispattive energy
        P_loadloss_pre=zeros(nbus-ngen,1);  % lossy dispattive energy
        Pe_mag_pre=zeros(ngen,1);   % Pe relevant to line inductance
        Qnet_mag_pre=zeros(nbus-ngen,1);   % Pe relevant to line inductance
        Pnet_mag_pre=zeros(nbus-ngen,1);   % Pe relevant to line inductance
        mag_load_P_pre=zeros(nbus-ngen,1);   % P Q load function on net
        mag_load_Q_pre=zeros(nbus-ngen,1);   % P Q load function on net
        for tm=2:no_duration
            ddelta_tmp=deltac_timedomain_DAE(n_start-1+tm,:)-deltac_timedomain_DAE(n_start-1+tm-1,:);
            dtheta_tmp=delta_net_clear_DAE(n_start-1+tm,:)-delta_net_clear_DAE(n_start-1+tm-1,:);
            dvoltage_tmp=voltage_net_clear_DAE(n_start-1+tm,:)-voltage_net_clear_DAE(n_start-1+tm-1,:);
            Pe=zeros(ngen,1);
            Pnet=zeros(nbus-ngen,1);
            Qnet=zeros(nbus-ngen,1);
            Pe_lossy=zeros(ngen,1);  % lossy dispattive energy
            Pnet_lossy=zeros(nbus-ngen,1);  % lossy dispattive energy
            Qnet_lossy=zeros(nbus-ngen,1);  % lossy dispattive energy
            P_loadloss=zeros(nbus-ngen,1);  % lossy dispattive energy
            Pe_mag=zeros(ngen,1);   % Pe relevant to line inductance
            Qnet_mag=zeros(nbus-ngen,1);   % Pe relevant to line inductance
            Pnet_mag=zeros(nbus-ngen,1);   % Pe relevant to line inductance

            mag_load_P=zeros(nbus-ngen,1);   % P Q load function on net
            mag_load_Q=zeros(nbus-ngen,1);   % P Q load function on net

            Ed_un=zeros(ngen,1);
            Ed_un_this=zeros(ngen,1);
            Ed_non=zeros(ngen,1);
            Ed_non_this=zeros(ngen,1);
            Pd_un=zeros(ngen,1);
            Pd_non=zeros(ngen,1);


            Ep_iter_gen=zeros(ngen,1);
            Ep_iter_net=zeros((nbus-ngen),1);


            % P calculation of gen
            for i=1:ngen
                for j=1:ngen
                    ddelta=deltac_timedomain_DAE(n_start-1+tm,i)-deltac_timedomain_DAE(n_start-1+tm,j);
                    Pe(i)=Pe(i)+E(i)*E(j)*B_post(i,j)*sin(ddelta)+E(i)*E(j)*G_post(i,j)*cos(ddelta);
                    if(i~=j)
                        Pe_lossy(i)=Pe_lossy(i)+E(i)*E(j)*G_post(i,j)*cos(ddelta);
                        Pe_mag(i)=Pe_mag(i)+E(i)*E(j)*B_post(i,j)*sin(ddelta);
                    end
                end
                for l=1:(nbus-ngen)
                    ddelta=deltac_timedomain_DAE(n_start-1+tm,i)-delta_net_clear_DAE(n_start-1+tm,l);
                    Pe(i)=Pe(i)+E(i)*voltage_net_clear_DAE(n_start-1+tm,l)*B_post(i,l+ngen)*sin(ddelta)+E(i)*voltage_net_clear_DAE(n_start-1+tm,l)*G_post(i,l+ngen)*cos(ddelta);
                    Pe_mag(i)=Pe_mag(i)+E(i)*voltage_net_clear_DAE(n_start-1+tm,l)*B_post(i,l+ngen)*sin(ddelta);
                    Pe_lossy(i)=Pe_lossy(i)+E(i)*voltage_net_clear_DAE(n_start-1+tm,l)*G_post(i,l+ngen)*cos(ddelta);
                end
                %damping energy
                Ed_un(i)=preset.d(i)*omegac_timedomain_DAE(n_start-1+tm-1,i)*omegab*ddelta_tmp(i);
                Ed_un_this(i)=preset.d(i)*omegac_timedomain_DAE(n_start-1+tm,i)*omegab*ddelta_tmp(i);
                Pd_un(i)=-preset.d(i)*(omegac_timedomain_DAE(n_start-1+tm-1,i)*omegab)^2;
                Ed_non(i)=preset.d(i)*(omegacoi_timedomain_DAE(n_start-1+tm-1,1)*omegab-omegab)*ddelta_tmp(i);
                Ed_non_this(i)=preset.d(i)*(omegacoi_timedomain_DAE(n_start-1+tm,1)*omegab-omegab)*ddelta_tmp(i);
                Pd_non(i)=-preset.d(i)*(omegacoi_timedomain_DAE(n_start-1+tm-1,1)*omegab-omegab)*(omegac_timedomain_DAE(n_start-1+tm-1,i)*omegab);
            end
            clear i j l
            % P calculation of Bus
            for i=1:(nbus-ngen)
                for j=1:ngen
                    ddelta=delta_net_clear_DAE(n_start-1+tm,i)-deltac_timedomain_DAE(n_start-1+tm,j);
                    Pnet(i)=Pnet(i)+voltage_net_clear_DAE(n_start-1+tm,i)*E(j)*B_post(i+ngen,j)*sin(ddelta)+voltage_net_clear_DAE(n_start-1+tm,i)*E(j)*G_post(i+ngen,j)*cos(ddelta);
                    Pnet_lossy(i)=Pnet_lossy(i)+voltage_net_clear_DAE(n_start-1+tm,i)*E(j)*G_post(i+ngen,j)*cos(ddelta);
                    Pnet_mag(i)=Pnet_mag(i)+voltage_net_clear_DAE(n_start-1+tm,i)*E(j)*B_post(i+ngen,j)*sin(ddelta);
                end
                for l=1:(nbus-ngen)
                    ddelta=delta_net_clear_DAE(n_start-1+tm,i)-delta_net_clear_DAE(n_start-1+tm,l);
                    Pnet(i)=Pnet(i)+voltage_net_clear_DAE(n_start-1+tm,i)*voltage_net_clear_DAE(n_start-1+tm,l)*B_post(i+ngen,l+ngen)*sin(ddelta)+voltage_net_clear_DAE(n_start-1+tm,i)*voltage_net_clear_DAE(n_start-1+tm,l)*G_post(i+ngen,l+ngen)*cos(ddelta);
                    if(i~=l)
                        Pnet_lossy(i)=Pnet_lossy(i)+voltage_net_clear_DAE(n_start-1+tm,i)*voltage_net_clear_DAE(n_start-1+tm,l)*G_post(i+ngen,l+ngen)*cos(ddelta);
                        Pnet_mag(i)=Pnet_mag(i)+voltage_net_clear_DAE(n_start-1+tm,i)*voltage_net_clear_DAE(n_start-1+tm,l)*B_post(i+ngen,l+ngen)*sin(ddelta);
                    else
                        P_loadloss(i) = voltage_net_clear_DAE(n_start-1+tm,i)^2*G_post(i+ngen,i+ngen);
                    end
                end
                for h=1:size(preset.Sload,1)
                    if (preset.Sload(h,1)==postfault.Transform(i+ngen))
                         Pnet(i)=Pnet(i)+preset.Sload(h,2);
                         mag_load_P(i)=mag_load_P(i)+preset.Sload(h,2);
                    end
                end
            end
            clear i j l h
            % Q/V calculation of Bus
            for i=1:(nbus-ngen)
                for j=1:ngen
                    ddelta=delta_net_clear_DAE(n_start-1+tm,i)-deltac_timedomain_DAE(n_start-1+tm,j);
                    Qnet(i)=Qnet(i)-E(j)*B_post(i+ngen,j)*cos(ddelta)+E(j)*G_post(i+ngen,j)*sin(ddelta);
                    Qnet_mag(i)=Qnet_mag(i)-E(j)*B_post(i+ngen,j)*cos(ddelta);
                    Qnet_lossy(i)=Qnet_lossy(i)+E(j)*G_post(i+ngen,j)*sin(ddelta);
                end
                for l=1:(nbus-ngen)
                    ddelta=delta_net_clear_DAE(n_start-1+tm,i)-delta_net_clear_DAE(n_start-1+tm,l);
                    Qnet(i)=Qnet(i)-voltage_net_clear_DAE(n_start-1+tm,l)*B_post(i+ngen,l+ngen)*cos(ddelta)+voltage_net_clear_DAE(n_start-1+tm,l)*G_post(i+ngen,l+ngen)*sin(ddelta);
                    Qnet_mag(i)=Qnet_mag(i)-voltage_net_clear_DAE(n_start-1+tm,l)*B_post(i+ngen,l+ngen)*cos(ddelta);
                    Qnet_lossy(i)=Qnet_lossy(i)+voltage_net_clear_DAE(n_start-1+tm,l)*G_post(i+ngen,l+ngen)*sin(ddelta);
                end
                for h=1:size(preset.Sload,1)
                    if (preset.Sload(h,1)==postfault.Transform(i+ngen))
                           Qnet(i)=Qnet(i)+preset.Sload(h,3)/voltage_net_clear_DAE(n_start-1+tm,i);
                           mag_load_Q(i)=mag_load_Q(i)+preset.Sload(h,3)/voltage_net_clear_DAE(n_start-1+tm,i);
                    end
                end
            end  
            clear i j l h

            if tm==2
                Pe_pre=Pe;
                Pnet_pre=Pnet;
                Qnet_pre=Qnet;
                Pe_lossy_pre=Pe_lossy;  % lossy dispattive energy
                Pnet_lossy_pre=Pnet_lossy;  % lossy dispattive energy
                Qnet_lossy_pre=Qnet_lossy;  % lossy dispattive energy
                P_loadloss_pre=P_loadloss;  % lossy dispattive energy
                Pe_mag_pre=Pe_mag;   % Pe relevant to line inductance
                Qnet_mag_pre=Qnet_mag;   % Pe relevant to line inductance
                Pnet_mag_pre=Pnet_mag;   % Pe relevant to line inductance
                mag_load_P_pre=mag_load_P;
                mag_load_Q_pre=mag_load_Q;
            end


            for i=1:ngen
                Ep_iter_gen(i)=(Pm(i)-Pe(i)+Pm(i)-Pe_pre(i))/2*ddelta_tmp(i);
                Ek(tm)=Ek(tm)+0.5*m(i)*(omegac_timedomain_DAE(n_start+tm-1,i)*omegab)^2;
            end
            for i=1:nbus-ngen
                Ep_iter_net(i)=(Pnet(i)+Pnet_pre(i))/2*dtheta_tmp(i) + (Qnet(i)+Qnet_pre(i))/2*dvoltage_tmp(i);
            end
            

            Ep_lossy_iter(tm)=Ep_lossy_iter(tm-1)+ddelta_tmp*(Pe_lossy+Pe_lossy_pre)/2+dtheta_tmp*(Pnet_lossy+Pnet_lossy_pre)/2+dvoltage_tmp*(Qnet_lossy+Qnet_lossy_pre)/2;   %losses energy on network exclude self-susceptance
            Ep_mag_iter(tm)=Ep_mag_iter(tm-1)+ddelta_tmp*(Pe_mag+Pe_mag_pre)/2+dtheta_tmp*(Pnet_mag+Pnet_mag_pre)/2+dvoltage_tmp*(Qnet_mag+Qnet_mag_pre)/2;% %magnatic potential energy
            Ep_loadloss_iter(tm)=Ep_loadloss_iter(tm-1)+dtheta_tmp*(P_loadloss+P_loadloss_pre)/2; %R-type load losses containing self-susceptance
            
            Ep_iter_record(tm)=Ep_iter_record(tm-1)-sum(Ep_iter_gen) + sum(Ep_iter_net); %sum Ep Iteration
            Ep_iter_record_withoutnet(tm)= Ep_iter_record(tm-1)-sum(Ep_iter_gen); %sum Ep Iteration

            Ep_magload_iter(tm) = Ep_magload_iter(tm-1)+dtheta_tmp*(mag_load_P+mag_load_P_pre)/2+dvoltage_tmp*(mag_load_Q_pre+mag_load_Q)/2;



            Ed_non_record(tm)=Ed_non_record(tm-1)+(sum(Ed_non)+sum(Ed_non_this))/2;
            Ed_un_record(tm)=Ed_un_record(tm-1)+(sum(Ed_un)+sum(Ed_un_this))/2;

            Err_step(tm)=Ep_iter_record(tm)+Ed_un_record(tm)+Ed_non_record(tm)+Ek(tm)-Ek(1); %sum energy
            Err_step_exp(tm)=Ep(tm)+Ed_un_record(tm)+Ed_non_record(tm)+Ek(tm)-Ek(1); %sum energy

            P_un_record(tm)=sum(Pd_un);
            P_non_record(tm)=sum(Pd_non);

            % Trapezoidal Int
            Pe_pre=Pe;
            Pnet_pre=Pnet;
            Qnet_pre=Qnet;
            Pe_lossy_pre=Pe_lossy;  % lossy dispattive energy
            Pnet_lossy_pre=Pnet_lossy;  % lossy dispattive energy
            Qnet_lossy_pre=Qnet_lossy;  % lossy dispattive energy
            P_loadloss_pre=P_loadloss;  % lossy dispattive energy
            Pe_mag_pre=Pe_mag;   % Pe relevant to line inductance
            Qnet_mag_pre=Qnet_mag;   % Pe relevant to line inductance
            Pnet_mag_pre=Pnet_mag;   % Pe relevant to line inductance
            mag_load_P_pre=mag_load_P;
            mag_load_Q_pre=mag_load_Q;
        end                                                  

    %% plot

    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end),Ep_loadloss,'color','r','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Ep_loadloss_iter,'color','g','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Potential Ray','Iter');
    title('Resistive lossy energy on load');
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end),Ep_lossy,'color','r','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Ep_lossy_iter,'color','g','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Calculation','Iter');
    title('lossy energy in grid');
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end),Ep_mag,'color','r','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Ep_mag_iter,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Potential','Iter');
    title('magnetic energy');
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end),Ep_mag_load,'color','r','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Ep_magload_iter,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Potential','Iter');
    title('magnetic load energy');
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end), Ep_loadloss-Ep_loadloss_iter,'color','g','LineWidth',2);  hold on; 
    plot(t_timedomain_DAE(n_start:n_end), Ep_lossy-Ep_lossy_iter,'color','y','LineWidth',2,'LineStyle','-');  hold on;
    plot(t_timedomain_DAE(n_start:n_end), Ep-Ep_iter_record,'color','b','LineWidth',2,'LineStyle','-');  hold on;
    plot(t_timedomain_DAE(n_start:n_end), Ep_loadloss-Ep_loadloss_iter+Ep_lossy-Ep_lossy_iter,'color','r','LineWidth',2,'LineStyle','--');  hold on;
    plot(t_timedomain_DAE(n_start:n_end), Ep_iter_record-Ep_iter_record_withoutnet,'color','k','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('load Lossy part err','network Lossy part err','Total err','Total lossy err','network potential error');
    title('Potential Energy Error');  

    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end),Ed_un_record,'color','r','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Ed_non_record,'color','b','LineWidth',2,'LineStyle','-');  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Ed_non_record+Ed_un_record,'color','k','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('uniformed damping','non-uniformed damping','total');
    title('Damping Energy');
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end),Ep+Ek,'color','r','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Ep_iter_record+Ek,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Expression','Iteration');
    title('Energy Function'); 
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end),P_un_record,'color','r','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),P_non_record,'color','g','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),P_non_record+P_un_record,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Uniform','Non-uniform','total');
    title('Derivative of damping energy'); 
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end),Err_step_exp,'color','r','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Err_step,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Expression','Iteration');
    title('Total Energy Change');  


%% ode function
function dfdt = f_timedomain(t,x)
    dfdt = F_3M9B_MR_ODE(x);
end

function dfdt = f_timedomain_DAE(t,x,system)
    dfdt = F_3M9B_SP_DAE(x,system);
end
function dfdt = f_timedomain_DAE_withflove(t,x,net_guess,system)
    net_value=fsolve(@(x)Fun_AEfslove_SPM(x,delta0,preset,system),net_guess,optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter','TolX',1e-9));
    dfdt = F_3M9B_SP_ODE(x,net_value,system);
end










