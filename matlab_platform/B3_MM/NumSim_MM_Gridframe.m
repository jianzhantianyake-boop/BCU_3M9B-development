% =========================================================================
% 教学操作说明：网络约简模型的故障前、故障中和故障后数值轨迹及能量绘图脚本。
% 使用方法：
%   新手只通过 run_bcu_beginner.m 的 "reduced_numerical" 模式运行；它会先执行
%   Cal_MM_CCT.m。手动运行本文件前，必须已在同一 base workspace 成功得到
%   prefault、fault、postfault 和 postfault.CUEP_delta。
% 参数：
%   无函数形参。Iter.Tfault（故障施加时刻）、Iter.Trecover（故障切除时刻）、
%   Iter.Ttotal（总时长）、Iter.Tunit（积分步长）均在下方设置且单位为 s；需满足
%   0 <= Tfault < Trecover <= Ttotal。PathEnergyCal 的多次赋值用于原始能量方法比较。
% 返回 / 工作区结果：
%   在 base workspace 生成/更新 IterData、能量比较变量及原始图窗，包括相角、
%   角速度、相平面和能量曲线。具体字段以运行完成后的 whos IterData 为准。
% 步骤：
%   1. 读取 CCT/初始化结果。2. 在故障前、故障中、故障后三段网络上推进原始
%   数值积分。3. 计算能量项和误差比较。4. 绘制原始时间域与相平面图。
% 单位：
%   Iter 时间变量为 s；转子角为 rad，角速度为 rad/s；功率、电压和导纳通常为 pu。
% 前置条件：
%   CCT 必须找到 CUEP；不要直接从空工作区运行，也不要以零向量伪造 CUEP。
%   减小 Iter.Tunit 或增大 Iter.Ttotal 会显著增加运行时间与内存占用。
% 研究与验证边界：
%   本次只更新注释，不改变积分、能量比较、时间设置或绘图逻辑；不同
%   PathEnergyCal 曲线的差异需要结合原始定义分析，不能直接视为精确误差证明。
% =========================================================================
%close all
%clear
%Cal_MM_Static;
%clear pfdata EMF Yload Case netdata

omegab=Basevalue.omegab;
%% 1) 时间设置：Tfault 是故障施加时刻，Trecover 是故障切除时刻。
    Iter.Tfault=20;
    Iter.Trecover=20.24;
    Iter.Ttotal=100;
    Iter.Tunit=1e-4;
    T_before = 10;
    T_after = 40;
%% 2) 固定步长轨迹：分别在 prefault/fault/postfault 三个 Yred 上推进 RK4。
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
%% 3) 自适应 ODE 对照：用同一三段工况生成时域参考轨迹。
% prefault
    system="prefault";
    delta0=prefault.SEP_delta;
    omega0=prefault.SEP_omegapu*omegab*ones(ngen,1);
    [t_prefault, x_prefault_all] = ode78(@f_timedomain,[0,Iter.Tfault],[delta0; omega0],odeset('RelTol',1e-5));
% faults
    system="fault";
    delta0=x_prefault_all(end,1:3)';
    omega0=x_prefault_all(end,4:6)';
    [t_fault, x_fault_all] = ode78(@f_timedomain,[Iter.Tfault,Iter.Trecover],[delta0; omega0],odeset('RelTol',1e-5));
% postfault
    system="postfault";
    delta0=x_fault_all(end,1:3)';
    omega0=x_fault_all(end,4:6)';
    [t_postfault, x_postfault_all] = ode78(@f_timedomain,[Iter.Trecover,Iter.Ttotal],[delta0; omega0],odeset('RelTol',1e-5));
    clear delta0 omega0
% data collection
t_timedomain = [t_prefault;t_fault;t_postfault];
deltac_timedomain = [x_prefault_all(:,1:3); x_fault_all(:,1:3); x_postfault_all(:,1:3)];
omega_timedomain = [x_prefault_all(:,4:6); x_fault_all(:,4:6); x_postfault_all(:,4:6)]./omegab;
omegacoi_timedomain = omega_timedomain*preset.m./sum(preset.m);
omegac_timedomain = omega_timedomain - omegacoi_timedomain*ones(1,3);


%% 4) 数据拼接：把三段轨迹拼成统一的 IterData，供能量和绘图使用。
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

    %% Trajectories in COI
    % thetac
        clear ylim
        figure;
        set(gca,'position',[0.115,0.12,0.815,0.84]);
        set(gcf,'position',[60 200 600 450]);

        n_start=cycle_pre-fix(T_before/Iter.Tunit)+1;
        n_end=cycle_pre+fix(T_after/Iter.Tunit);
        
        plot(IterData.Tout(n_start:n_end),IterData.thetac((n_start:n_end),1),'linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
        plot(IterData.Tout(n_start:n_end),IterData.thetac((n_start:n_end),2),'linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
        plot(IterData.Tout(n_start:n_end),IterData.thetac((n_start:n_end),3),'linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
        plot(IterData.Tout(n_start:n_end),IterData.deltacd((n_start:n_end),1),'linewidth',2,'color',[200/255 200/255 200/255]);    hold on;

        %ode time domain
        plot(t_timedomain,deltac_timedomain(:,1),'LineStyle',':','linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
        plot(t_timedomain,deltac_timedomain(:,2),'LineStyle',':','linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
        plot(t_timedomain,deltac_timedomain(:,3),'LineStyle',':','linewidth',2,'color',[0/255 175/255 0/255]);    hold on;

        xlim([20-T_before 20+T_after]);
        grid on;    grid minor;
        ylim([-2,4]);
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
        xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
        ylabel('\delta_c(rad)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);


        figure;
        set(gca,'position',[0.115,0.12,0.815,0.84]);
        set(gcf,'position',[60 200 600 450]);
        plot(IterData.Tout(n_start:n_end),IterData.omegac((n_start:n_end),1),'linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
        plot(IterData.Tout(n_start:n_end),IterData.omegac((n_start:n_end),2),'linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
        plot(IterData.Tout(n_start:n_end),IterData.omegac((n_start:n_end),3),'linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
        
        plot(IterData.Tout(n_start:n_end),IterData.omegacoi((n_start:n_end),1)-1,'linewidth',2,'color',[0/255 0/255 0/255]);    hold on;
        plot(IterData.Tout(n_start:n_end),IterData.omegacd((n_start:n_end),1),'linewidth',2,'color',[200/255 200/255 200/255]);    hold on;

        %ode time domain
        plot(t_timedomain,omegac_timedomain(:,1),'LineStyle',':','linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
        plot(t_timedomain,omegac_timedomain(:,2),'LineStyle',':','linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
        plot(t_timedomain,omegac_timedomain(:,3),'LineStyle',':','linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
        plot(t_timedomain,omegacoi_timedomain-1,'LineStyle',':','linewidth',2,'color',[0/255 0/255 0/255]);    hold on;

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
% 5) 能量后处理：由轨迹计算动能、势能、阻尼能和临界轨迹差异。
    no_duration = n_end - n_start + 1;
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    Yred_post=postfault.Yred;
    G_post=real(Yred_post);
    B_post=imag(Yred_post);
    % Potential energy
        [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3)]=Fun_Cal_PotentialEnergy(preset,postfault,postfault.SEP_delta,IterData.thetac(n_start,:)');
        Ep_start=sum(Ep_tmp);
        Ep_lossy_start=Ep_tmp(3);
        Ep_mag_start=Ep_tmp(2);
        clear Ep_tmp
        [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3)]=Fun_Cal_PotentialEnergy(preset,postfault,postfault.SEP_delta,IterData.thetac(cycle_pre+cycle_fault,:)');
        Ep_end=sum(Ep_tmp);   clear Ep_tmp
        [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3)]=Fun_Cal_PotentialEnergy(preset,postfault,postfault.SEP_delta,postfault.CUEP_delta);
        Ep_CUEP=sum(Ep_tmp);   clear Ep_tmp
        Ep=zeros(no_duration,1);
        Ep_path=zeros(no_duration,1);
        Ep_mag=zeros(no_duration,1);
        EP_trapezoidal=zeros(no_duration,1);
        EP_trapezoidal_multi=zeros(no_duration,1);
        preset.PathEnergyCal=0;
        for tm=1:size(Ep,1)
            [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3)]=Fun_Cal_PotentialEnergy(preset,postfault,postfault.SEP_delta,IterData.thetac(n_start-1+tm,:)');
            Ep_path(tm)=Ep_tmp(3);
            Ep_mag(tm)=Ep_tmp(2);
            Ep(tm)=sum(Ep_tmp);   clear Ep_tmp
        end
        preset.PathEnergyCal=1;
        for tm=1:size(Ep,1)
            [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3)]=Fun_Cal_PotentialEnergy(preset,postfault,postfault.SEP_delta,IterData.thetac(n_start-1+tm,:)');
            EP_trapezoidal(tm)=Ep_tmp(3);   clear Ep_tmp
        end
        preset.PathEnergyCal=10;
        for tm=1:size(Ep,1)
            [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3)]=Fun_Cal_PotentialEnergy(preset,postfault,postfault.SEP_delta,IterData.thetac(n_start-1+tm,:)');
            EP_trapezoidal_multi(tm)=Ep_tmp(3);   clear Ep_tmp
        end
        preset.PathEnergyCal=0;
     % Kinetic energy
        Ek0=zeros(ngen,1);
        Eke=zeros(ngen,1);
        for i=1:ngen
            Ek0(i)=0.5*m(i)*(IterData.omegac(n_start,i)*omegab)^2;
            Eke(i)=0.5*m(i)*(IterData.omegac(cycle_pre+cycle_fault,i)*omegab)^2;
        end
        Ek_start=sum(Ek0);
        Ek_end=sum(Eke);
    % Damping energy in iteration
        % for Potential energy error display
        Ed_non_record=zeros(no_duration,1);
        Ed_un_record=zeros(no_duration,1);
        Ep_lossy_iter=zeros(no_duration,1);
        Ep_mag_iter=zeros(no_duration,1);
        Ep_iter_record=zeros(no_duration,1);
        P_un_record=zeros(no_duration,1);
        P_non_record=zeros(no_duration,1);
        Ep_iter_record(1)=Ep_start;   
        Ep_lossy_iter(1)=Ep_lossy_start;
        Ep_mag_iter(1)=Ep_mag_start;
        Ek=zeros(no_duration,1);
        Ek(1)=Ek_start;
        Err_step=zeros(no_duration,1);
        Err_step_exp=zeros(no_duration,1);
        D_un_record=zeros(no_duration,1);

        Pe_pre=zeros(ngen,1);
        D_un_critical = preset.d'*(IterData.omegacoi(cycle_pre+cycle_fault+1,1)*omegab-omegab+postfault.SEP_omegapu*omegab-omegab)*(postfault.CUEP_delta-postfault.SEP_delta)/2;
        for tm=2:no_duration
            ddelta_tmp=IterData.thetac(n_start-1+tm,:)-IterData.thetac(n_start-1+tm-1,:);
            Pe=zeros(ngen,1);
            Pe_lossy=zeros(ngen,1);  % lossy dispattive energy
            Pe_mag=zeros(ngen,1);   % Pe relevant to line inductance
            Ed_un=zeros(ngen,1);
            Ed_un_this=zeros(ngen,1);
            Ed_non=zeros(ngen,1);
            Ed_non_this=zeros(ngen,1);
            Pd_un=zeros(ngen,1);
            Pd_non=zeros(ngen,1);
            Ep_iter_gen=zeros(ngen,1);
            D_non_gen=zeros(ngen,1);
            for i=1:ngen
                for j=1:ngen
                    delta=IterData.thetac(n_start-1+tm,i)-IterData.thetac(n_start-1+tm,j);
                    Pe(i)=Pe(i)+E(i)*E(j)*B_post(i,j)*sin(delta)+E(i)*E(j)*G_post(i,j)*cos(delta);
                    if(i~=j)
                    Pe_lossy(i)=Pe_lossy(i)+E(i)*E(j)*G_post(i,j)*cos(delta);
                    Pe_mag(i)=Pe_mag(i)+E(i)*E(j)*B_post(i,j)*sin(delta);
                    end
                end
            end

            if tm==2
                Pe_pre=Pe;
            end

            for i=1:ngen
                Ed_un(i)=preset.d(i)*IterData.omegac(n_start-1+tm-1,i)*omegab*ddelta_tmp(i);
                Ed_un_this(i)=preset.d(i)*IterData.omegac(n_start-1+tm,i)*omegab*ddelta_tmp(i);
                Pd_un(i)=-preset.d(i)*(IterData.omegac(n_start-1+tm-1,i)*omegab)^2;
                Ed_non(i)=preset.d(i)*(IterData.omegacoi(n_start-1+tm-1,1)*omegab-omegab)*ddelta_tmp(i);
                Ed_non_this(i)=preset.d(i)*(IterData.omegacoi(n_start-1+tm,1)*omegab-omegab)*ddelta_tmp(i);
                Pd_non(i)=-preset.d(i)*(IterData.omegacoi(n_start-1+tm-1,1)*omegab-omegab)*(IterData.omegac(n_start-1+tm-1,i)*omegab);
                Ep_iter_gen(i)=(Pm(i)-Pe(i)+Pm(i)-Pe_pre(i))/2*ddelta_tmp(i);
                Ek(tm)=Ek(tm)+0.5*m(i)*(IterData.omegac(n_start+tm-1,i)*omegab)^2;
                D_non_gen(i) = preset.d(i)*(IterData.omegacoi(n_start-1+tm-1,1)*omegab-omegab+postfault.SEP_omegapu*omegab-omegab)*(IterData.thetac(n_start-1+tm-1,i)-postfault.SEP_delta(i))/2;
            end
            Pe_pre = Pe;

            Ep_lossy_iter(tm)=Ep_lossy_iter(tm-1)+ddelta_tmp*Pe_lossy;   
            Ep_mag_iter(tm)=Ep_mag_iter(tm-1)+ddelta_tmp*Pe_mag;          
            Ep_iter_record(tm)=Ep_iter_record(tm-1)-sum(Ep_iter_gen);%-ddelta_tmp*(Pm-Pe);
            Ed_non_record(tm)=Ed_non_record(tm-1)+(sum(Ed_non)+sum(Ed_non_this))/2;
            Ed_un_record(tm)=Ed_un_record(tm-1)+(sum(Ed_un)+sum(Ed_un_this))/2;
            Err_step(tm)=Ep_iter_record(tm)+Ed_un_record(tm)+Ed_non_record(tm)+Ek(tm)-Ek(1);
            Err_step_exp(tm)=Ep(tm)+Ed_un_record(tm)+Ed_non_record(tm)+Ek(tm)-Ek(1);
            P_un_record(tm)=sum(Pd_un);
            P_non_record(tm)=sum(Pd_non);
            D_un_record(tm)=sum(D_non_gen);
        end                                                  

    %% plot

    figure;
    plot(IterData.Tout(n_start:n_end),Ep_path,'color','r','LineWidth',2);  hold on;
    plot(IterData.Tout(n_start:n_end),EP_trapezoidal,'color','g','LineWidth',2,'LineStyle','--');  hold on;
    plot(IterData.Tout(n_start:n_end),EP_trapezoidal_multi,'color','y','LineWidth',2,'LineStyle','--');  hold on;
    plot(IterData.Tout(n_start:n_end),Ep_lossy_iter,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Potential Ray','Potential trapezoidal','Potential multistep trapezoidal','Iter');
    title('lossy energy');
    figure;
    plot(IterData.Tout(n_start:n_end),Ep_mag,'color','r','LineWidth',2);  hold on;
    plot(IterData.Tout(n_start:n_end),Ep_mag_iter,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Potential','Iter');
    title('magnetic energy');
    figure;
    plot(IterData.Tout(n_start:n_end),Ep_path-Ep_lossy_iter,'color','r','LineWidth',2);  hold on;
    plot(IterData.Tout(n_start:n_end),Ep-Ep_iter_record,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Lossy part err','Total err');
    title('Error');  
    figure;
    plot(IterData.Tout(n_start:n_end),Ed_un_record,'color','r','LineWidth',2);  hold on;
    plot(IterData.Tout(n_start:n_end),Ed_non_record,'color','b','LineWidth',2,'LineStyle','-');  hold on;
    plot(IterData.Tout(n_start:n_end),Ed_non_record+Ed_un_record,'color','k','LineWidth',2,'LineStyle','--');  hold on;
    plot(IterData.Tout(n_start:n_end),D_un_record,'color','m','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('uniformed damping','non-uniformed damping','total','non-uniformed damping energy estimation');
    title('Damping Energy');
    figure;
    plot(IterData.Tout(n_start:n_end),Ep+Ek,'color','r','LineWidth',2);  hold on;
    plot(IterData.Tout(n_start:n_end),Ep_iter_record+Ek,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    plot(IterData.Tout(n_start:n_end),Ep_CUEP*ones(size(IterData.Tout(n_start:n_end))),'color','k','LineWidth',2,'LineStyle','--');  hold on;
    plot(IterData.Tout(n_start:n_end),Ep+Ek+D_un_record,'color','m','LineWidth',2,'LineStyle','--');  hold on;
    plot(IterData.Tout(n_start:n_end),(Ep_CUEP+D_un_critical)*ones(size(IterData.Tout(n_start:n_end))),'color','k','LineWidth',2,'LineStyle',':');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Expression','Iteration','Critical energy','Revised Energy Function','Revised Critical energy');
    title('Energy Function'); 
    figure;
    plot(IterData.Tout(n_start:n_end),P_un_record,'color','r','LineWidth',2);  hold on;
    plot(IterData.Tout(n_start:n_end),P_non_record,'color','g','LineWidth',2);  hold on;
    plot(IterData.Tout(n_start:n_end),P_non_record+P_un_record,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Uniform','Non-uniform','total');
    title('Derivative of damping energy'); 
    figure;
    plot(IterData.Tout(n_start:n_end),Err_step_exp,'color','r','LineWidth',2);  hold on;
    plot(IterData.Tout(n_start:n_end),Err_step,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Expression','Iteration');
    title('Total Energy Change');  
    figure;
    plot(IterData.deltacd((n_start:cycle_pre+cycle_fault),1),IterData.omegacoi((n_start:cycle_pre+cycle_fault),1)-1,'color','r','LineWidth',2);  hold on;
    plot(IterData.deltacd((cycle_pre+cycle_fault+1:n_end),1),IterData.omegacoi((cycle_pre+cycle_fault+1:n_end),1)-1,'color','b','LineWidth',2);  hold on;
    xlabel('\delta_D');
    ylabel('\omega_{COI}');
    legend('fault-on','post-fault');
    title('Damping process'); 
  
%% ode function
function dfdt = f_timedomain(t,x)
    dfdt = F_3M9B_MR_ODE(x);
end










