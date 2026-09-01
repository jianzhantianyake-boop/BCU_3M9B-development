% =========================================================================
% 操作说明：网络约简模型的 CCT、逃逸点、MGP 与 CUEP 计算脚本。
% 使用方法：
%   从 run_bcu_beginner.m 选择 "reduced_cct" 或 "reduced_numerical"。
%   操作时，在项目根目录执行 run(fullfile('B3_MM','Cal_MM_CCT.m'))；
%   本脚本会 clear 后自行调用 Cal_MM_Static.m，因此不要预先塞入同名工作区变量。
% 参数：
%   无函数形参。此文件顶部的 Tfault 是故障轨迹积分时长（s），Tunit 是积分
%   步长（s）；机器、故障和 case 参数应在 Cal_MM_Static.m 的设置段修改。
% 返回 / 工作区结果：
%   在 base workspace 生成或更新 escape、fault.traj、MGP、postfault.CUEP_delta
%   等结果，并按原实现绘制故障轨迹、SEP、MGP、CUEP 的相平面图。
% 步骤：
%   1. 清理工作区并完成网络约简初始化。2. 沿故障网络积分得到逃逸点。
%   3. 搜索 MGP/CUEP 并计算能量相关判据。4. 以原始图形标记关键状态点。
% 单位：
%   Tfault、Tunit 为 s；相角为 rad，角速度为 rad/s；网络、电势与功率通常为 pu。
% 前置条件：
%   需要路径、MATPOWER、fsolve 和原始数值函数可用。脚本会 clear，随后生成的
%   CUEP 是 NumSim_MM_Gridframe.m 的前置数据；未找到 CUEP 时不要继续跑轨迹。
% =========================================================================
%close all
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cal_MM_Static;
clear EMF Yload Case netdata
clear DT HT ngen
%% 1) 故障轨迹退出点：从预故障 SEP 沿故障网络积分，记录边界退出状态。
    BCUCFG = bcu_config();
    Tfault=bcu_pick(BCUCFG,'Tfault',2);   Tunit=bcu_pick(BCUCFG,'Tunit',1e-4);
    delta0=prefault.SEP_delta;
    omega0=prefault.SEP_omegapu*Basevalue.omegab;
    [theta_RK4,omega_RK4,thetac_RK4,omegac_RK4,escape.tm]=Fun_Cal_Exitpoint(Tfault,Tunit,fault.Yred,postfault.Yred,preset,delta0,omega0,Basevalue.omegab);
    escape.thetac=thetac_RK4(escape.tm,:);
    escape.omegac=omegac_RK4(escape.tm,:);
    escape.theta=theta_RK4(escape.tm,:);
    escape.omega=omega_RK4(escape.tm,:)-Basevalue.omegab;
    fault.traj.theta=theta_RK4;
    fault.traj.omega=omega_RK4;   % unit: rad/s    
    fault.traj.thetac=thetac_RK4;
    fault.traj.omegac=omegac_RK4;  % unit: rad/s
    fault.traj.Tunit=Tunit;
    fault.traj.Tlength=Tfault;
%     % exit point overwrite by critical value (get from numiter)
%     escape.thetac=[0.129933252577263	0.679495354275216	0.564006744715918	0.576040301569986	0.665736000959453	0.598079609167577	0.638117660815637	0.793490807600584	0.950721163419743	-0.339524838721402];
%     escape.theta=[1.86491858030534	2.41448068200330	2.29899207244400	2.31102562929807	2.40072132868753	2.33306493689566	2.37310298854372	2.52847613532867	2.68570649114782	1.39546048900668];
%     escape.omega=[343.145203761379	345.100261615261	344.755690005566	345.844936412527	343.957111976891	345.344513164961	345.651241344674	350.831688442615	348.797293319823	337.936771756932]-Basevalue.omegab;
%     ngen=size(prefault.Yred,1);
%     escape.omegac=escape.omega+(-340.763031348808+Basevalue.omegab)*ones(1,ngen);

    strEXIT=['Exit point is [' repmat('%1.4f ',1,numel(escape.thetac)) '] (in COI frame)\n'];
    fprintf(strEXIT,escape.thetac);
    f1=figure(1);
    f2=figure(2);
    figure(f2);
    grid on; hold on;
    %close(f1);
    %clear f1
    figure(f1);
    xlabel('\delta_2');
    ylabel('\delta_3');
    plot(escape.thetac(2),escape.thetac(3),'xr','LineWidth',1.5,'MarkerSize',10);
    grid on;
    hold on;
    plot(prefault.SEP_delta(2),prefault.SEP_delta(3),'.k','LineWidth',2,'MarkerSize',10);
    plot(fault.traj.thetac(:,2),fault.traj.thetac(:,3),'-','LineWidth',1.5,'color',[200/255 200/255 200/255]);
    axis([0,2.5,0,3.5]);
    %clear omega_RK4 omegac_RK4 theta_RK4 thetac_RK4
    clear Tfault delta0 omega0
%% 2) MGP：从退出点集合追踪故障后稳定边界，选择势场范数较小的点。
    [MGP.thetac_MGP,MGP.num_Traj,MGP.flag_MGP,Normtt, norm_min]=Fun_Cal_MGP(escape.thetac,postfault,preset);
    strMGP=['Selected MGP is [' repmat('%1.4f ',1,numel(MGP.thetac_MGP)) ']\n'];
    figure(f1);
    plot(MGP.thetac_MGP(2),MGP.thetac_MGP(3),'xr','LineWidth',1.5,'MarkerSize',10);
    fprintf(strMGP,MGP.thetac_MGP);
%% 3) CUEP：把 MGP 作为 fsolve 初值，寻找满足故障后平衡方程的临界点。
    %% use Newton powerflow method================
    if(preset.EquCal==1)
        [postfault.CUEP_delta,postfault.CUEP_omegapu,flag_iter,n_iter]=Fun_SEPiteration(postfault.Yred,preset.Pmpu,preset.Epu,preset.m,preset.d,MGP.thetac_MGP',0,Basevalue.omegab,1e5,1e-5);
    %% use fsolve function================
    elseif(preset.EquCal==2)
        ngen=size(prefault.Yred,1);
        deltaomega_init=zeros(ngen,1);
        for i=1:ngen-1
            deltaomega_init(i)=MGP.thetac_MGP(i)-MGP.thetac_MGP(ngen);
        end
        deltaomega_init(ngen)=0;    %from lei (?reason)
        Results_fsolve=fsolve(@(delta_omega)Fun_SEPfslove(delta_omega,preset,postfault,Basevalue),deltaomega_init,options);
        postfault.CUEP_omegapu=Results_fsolve(ngen)/Basevalue.omegab+1;
        delta_tmp=[Results_fsolve(1:ngen-1);0];
        postfault.CUEP_delta=delta_tmp-delta_tmp'*preset.m/sum(preset.m);
        clear delta_tmp Results_fsolve deltaomega_init
    end
        [postfault.CUEP_Perr,flag_SEPerr]=Fun_SEPcheck(postfault,preset,postfault.CUEP_delta,(postfault.CUEP_omegapu-1)*Basevalue.omegab);
    
%         for i=1:ngen-1
%             deltaomega_init(i)=postfault.CUEP_delta(i)-postfault.CUEP_delta(ngen);
%         end
%         deltaomega_init(ngen)=(postfault.CUEP_omegapu-1)*Basevalue.omegab;
%     Results_fsolve=fsolve(@(delta_omega)Fun_SEPfslove(delta_omega,preset,postfault,Basevalue),deltaomega_init,options);
%     omega_tmp=Results_fsolve(ngen)/Basevalue.omegab+1;
%     delta_tmp=[Results_fsolve(1:ngen-1)';0];
%     thetac_tmp=delta_tmp-delta_tmp'*preset.m/sum(preset.m);
%     [postfault.CUEP_Perr,flag_SEPerr]=Fun_SEPcheck(postfault,preset,thetac_tmp,(omega_tmp-1)*Basevalue.omegab);

    if(flag_SEPerr==1)
        error('SEP calculation error: the SEP solution cannot satisfy equilibrium condition!')
    end
    % check whether CUEP calculation is correct
    if(norm(postfault.CUEP_delta-postfault.SEP_delta)<1e-2)
        error('CUEP calculation error!');
    end
    strCUEP=['CUEP is [' repmat('%1.4f ',1,numel(postfault.CUEP_delta)) '] (in COI frame)\n'];
    fprintf('CUEP found!\n');
    fprintf(strCUEP, postfault.CUEP_delta);
    clear flag_iter n_iter
    clear strCUEP strMGP strEXIT
%     close stability following figure
    figure(f1);
    plot(postfault.CUEP_delta(2),postfault.CUEP_delta(3),'om','LineWidth',1.5,'MarkerSize',8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4) CCT：用 CUEP 的临界能量与故障轨迹能量交点估计切除时刻。
%     Critical Energy at CUEP
    [Ep(1),Ep(2),Ep(3)]=Fun_Cal_PotentialEnergy(preset,postfault,postfault.SEP_delta,postfault.CUEP_delta);
    E_critical=sum(Ep);
    % calculate CCT based on direct method
    [Critical.LEA.CCT,Critical.LEA.Exit_thetac,Critical.LEA.Exit_omegac,Critical.LEA.Exit_theta,Critical.LEA.Exit_omega,Critical.LEA.flag_CCT]=Fun_Cal_CCT_Energy(E_critical,fault,postfault,preset);
    Critical.Ep=E_critical;
    clear E_critical Ep
    strLyaCCT=['CCT(LEA) is ' repmat('%1.4f ',1,numel(Critical.LEA.CCT)) 's \n'];
    fprintf(strLyaCCT,Critical.LEA.CCT);
    clear strLyaCCT

    %% Calculate Real CCT
    [Critical.REA.CCT,Critical.REA.Exit_thetac,Critical.REA.Exit_omegac,Critical.REA.Exit_theta,Critical.REA.Exit_omega,Critical.REA.flag_CCT,Critical.Traj.Stb,Critical.Traj.Unstb]=Fun_Cal_CCT_Real(fault,postfault,preset,Basevalue,Critical.LEA.CCT);
    strREACCT=['CCT(REA) is ' repmat('%1.4f ',1,numel(Critical.REA.CCT)) 's \n'];
    fprintf(strREACCT,Critical.REA.CCT);
    clear strREACCT

    %% Estimate damping energy
%     Damping=zeros(ngen,1);
%     Critical.LEA.Exit_omegacoi = Critical.LEA.Exit_omega*preset.m/sum(preset.m)/Basevalue.omegab - 1;
%     for i=1:ngen
%         Damping(i)=preset.d(i)*((postfault.CUEP_omegapu-1)*Basevalue.omegab + Critical.LEA.Exit_omegacoi*Basevalue.omegab )*(postfault.CUEP_delta(i)-Critical.LEA.Exit_thetac(i));
%     end
%     Damping_est = sum(Damping);
%     E_critical=Critical.Ep + Damping_est;
% 
% 
%     [Critical.DEA2.CCT,Critical.DEA2.Exit_thetac,Critical.DEA2.Exit_omegac,Critical.DEA2.Exit_theta,Critical.DEA2.Exit_omega,Critical.DEA2.flag_CCT]=Fun_Cal_CCT_Energy(E_critical,fault,postfault,preset);
%     
%     clear E_critical Ep
%     strLyaCCT=['CCT(DEA2) is ' repmat('%1.4f ',1,numel(Critical.LEA.CCT)) 's \n'];
%     fprintf(strLyaCCT,Critical.DEA2.CCT);
%     clear strLyaCCT


%% Critical Stable and Unstable trajectories
%    [Energy,Group]=Fun_Cal_DampingEnergy(postfault,preset,Basevalue,Critical);
%    E_critical=Critical.Ep+Energy.Ed_uniform+Energy.Ed_nonuniform;
%    [Critical.DEA.CCT,Critical.DEA.Exit_thetac,Critical.DEA.Exit_omegac,Critical.DEA.Exit_theta,Critical.DEA.Exit_omega,Critical.DEA.flag_CCT]=Fun_Cal_CCT_Energy(E_critical,fault,postfault,preset);
%    
%    strDEACCT=['CCT(DEA) is ' repmat('%1.4f ',1,numel(Critical.DEA.CCT)) 's \n'];
%    fprintf(strDEACCT,Critical.DEA.CCT);
%    clear strDEACCT
% 
%    for i=1:20e4
%        dis_CUEP(i)=norm(Critical.Traj.Stb.thetac(i,:)-postfault.CUEP_delta');
%    end
%    figure(12);
%    plot(0:Tunit:Tunit*(20e4-1),dis_CUEP);
