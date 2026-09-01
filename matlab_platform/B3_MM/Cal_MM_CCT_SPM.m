% =========================================================================
% 操作说明：结构保持模型（SPM）的 CCT、逃逸点、MGP 与 CUEP 计算脚本。
% 使用方法：
%   从 run_bcu_beginner.m 选择 "spm_cct" 或 "spm_numerical"。操作时，
%   在项目根目录执行 run(fullfile('B3_MM','Cal_MM_CCT_SPM.m'))；本脚本会 clear
%   后自行调用 Cal_MM_Static_SPM.m，因此不要依赖运行前残留的 base workspace。
% 参数：
%   无函数形参。顶部 Tfault 为故障轨迹积分时长（s），Tunit 为积分步长（s）；
%   机器、ZIP 负荷、case 和故障配置应在 Cal_MM_Static_SPM.m 的设置段修改。
% 返回 / 工作区结果：
%   在 base workspace 生成/更新 escape、fault.traj、MGP、postfault.CUEP_delta
%   及 SPM 网络状态，并按原实现绘制关键相平面图。
% 步骤：
%   1. 清理工作区并完成 SPM 初始化。2. 沿故障状态积分确定逃逸点。3. 搜索
%   MGP/CUEP 与能量临界量。4. 输出原始状态点和相平面标记。
% 单位：
%   Tfault、Tunit 为 s；相角为 rad，角速度为 rad/s；电压、功率、导纳通常为 pu。
% 前置条件：
%   需要路径、MATPOWER、fsolve 和结构保持模型的原始数值函数。若没有得到
%   postfault.CUEP_delta，不能直接继续运行 NumSim_MM_Gridframe_SPM.m。
% =========================================================================
close all
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cal_MM_Static_SPM;
clear EMF Yload Case netdata
clear DT HT
%% 1) 结构保持故障轨迹退出点：状态同时包含网络角度和电压。
    Tfault=0.5;   Tunit=1e-4;
    delta0=prefault.SEP_delta;
    omega0=prefault.SEP_omegapu*Basevalue.omegab;
    delta_net0=prefault.net_delta;
    voltage_net0=prefault.net_voltage;
    [deltac,omega,omegac,theta,voltage,escape.tm,escape.Dotproduct]=Fun_Cal_Exitpoint_SPM(Tfault,Tunit,fault,postfault,preset,delta0,omega0,delta_net0,voltage_net0,Basevalue);
    escape.deltac=deltac(escape.tm,:);
    escape.omegac=omegac(escape.tm,:);
    escape.omega=omega(escape.tm,:)-Basevalue.omegab;
    escape.theta=theta(escape.tm,:);
    escape.voltage=voltage(escape.tm,:);

    fault.traj.deltac=deltac;
    fault.traj.omega=omega;   % unit: rad/s    
    fault.traj.theta=theta;
    fault.traj.voltage=voltage;
    fault.traj.omegac=omegac;  % unit: rad/s
    fault.traj.Tunit=Tunit;
    fault.traj.Tlength=Tfault;


    strEXIT=['Exit point is [' repmat('%1.4f ',1,numel(escape.deltac)) '] (in COI frame)\n'];
    fprintf(strEXIT,escape.deltac);
    f1=figure(1);
    f2=figure(2);
    figure(f2);
    grid on; hold on;
    figure(f1);
    xlabel('\delta_2');
    ylabel('\delta_3');
    plot(escape.deltac(2),escape.deltac(3),'xr','LineWidth',1.5,'MarkerSize',10);
    grid on;
    hold on;
    plot(prefault.SEP_delta(2),prefault.SEP_delta(3),'.k','LineWidth',2,'MarkerSize',10);
    plot(fault.traj.deltac(:,2),fault.traj.deltac(:,3),'-','LineWidth',1.5,'color',[200/255 200/255 200/255]);
    axis([0,2.5,0,3.5]);
    clear omega_RK4 omegac_RK4 theta_RK4 thetac_RK4
    clear Tfault delta0 omega0
%% 2) 结构保持 MGP：沿保留网络节点的边界轨迹寻找候选点。
    [MGP.detac_MGP,MGP.theta_MGP,MGP.voltage_MGP,MGP.num_Traj,MGP.flag_MGP,Normtt, norm_min]=Fun_Cal_MGP_SPM(escape.deltac,escape.theta,escape.voltage,postfault,preset);
    strMGP=['Selected MGP is [' repmat('%1.4f ',1,numel(MGP.detac_MGP)) ']\n'];
    figure(f1);
    plot(MGP.detac_MGP(2),MGP.detac_MGP(3),'xr','LineWidth',1.5,'MarkerSize',10);
    fprintf(strMGP,MGP.detac_MGP);
%% 3) 结构保持 CUEP：fsolve 同时满足发电机和网络节点平衡方程。
    %% use Newton powerflow method================
    if(preset.EquCal==1)
        %
    %% use fsolve function================
    elseif(preset.EquCal==2)
        ngen=size(prefault.Yred,1);
        deltaomega_init=zeros(ngen,1);
        thetanet_init=zeros((nbus-ngen),1);
        for i=1:ngen-1
            deltaomega_init(i)=MGP.detac_MGP(i)-MGP.detac_MGP(ngen);
        end
        for i=1:(nbus-ngen)
            thetanet_init(i)=MGP.theta_MGP(i)-MGP.detac_MGP(ngen);
        end
        deltaomega_init(ngen)=0;
        x_init=[deltaomega_init; thetanet_init; MGP.voltage_MGP'];
        options = optimset('TolFun',1e-50,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter','TolX',1e-9);
        options.StepTolerance = 1e-50;
        Results_fsolve=fsolve(@(x)Fun_SEPfslove_SPM(x,preset,postfault,Basevalue),x_init,options);
        postfault.CUEP_omegapu=Results_fsolve(ngen)/Basevalue.omegab+1;
        delta_tmp=[Results_fsolve(1:ngen-1);0];
        postfault.deltacoi=delta_tmp'*preset.m/sum(preset.m);
        postfault.CUEP_delta=delta_tmp-postfault.deltacoi;
        postfault.CUEP_net_theta=Results_fsolve((ngen+1):nbus)-deltacoi;
        postfault.CUEP_net_voltage=Results_fsolve((nbus+1):(2*nbus-ngen));
        %clear delta_tmp Results_fsolve deltaomega_init
    end
        [postfault.CUEP_Perr,flag_SEPerr]=Fun_SEPcheck(postfault,preset,postfault.CUEP_delta,(postfault.CUEP_omegapu-1)*Basevalue.omegab);
    
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
%% 4) 结构保持能量 CCT：用扩展状态轨迹计算临界能量交点。
%     Critical Energy at CUEP
    clear Ep
    preset.PathEnergyCal=0;
    [Ep(1),Ep(2),Ep(3),Ep(4),Ep(5)]=Fun_Cal_PotentialEnergy_SPM(preset,postfault,postfault.CUEP_delta,postfault.CUEP_net_theta,postfault.CUEP_net_voltage);
    E_critical=sum(Ep);
    % calculate CCT based on direct method
    [Critical.LEA.CCT,Critical.LEA.Exit_detac,Critical.LEA.Exit_omegac,Critical.LEA.Exit_theta,Critical.LEA.Exit_omega,Critical.LEA.Exit_voltage,Critical.LEA.flag_CCT]=Fun_Cal_CCT_Energy_SPM(E_critical,fault,postfault,preset);
    Critical.Ep=E_critical;
    clear E_critical Ep
    strLyaCCT=['CCT(LEA) is ' repmat('%1.4f ',1,numel(Critical.LEA.CCT)) 's \n'];
    fprintf(strLyaCCT,Critical.LEA.CCT);
    clear strLyaCCT

%     %% Calculate Real CCT
%     [Critical.REA.CCT,Critical.REA.Exit_thetac,Critical.REA.Exit_omegac,Critical.REA.Exit_theta,Critical.REA.Exit_omega,Critical.REA.flag_CCT,Critical.Traj.Stb,Critical.Traj.Unstb]=Fun_Cal_CCT_Real(fault,postfault,preset,Basevalue,Critical.LEA.CCT);
%     strREACCT=['CCT(REA) is ' repmat('%1.4f ',1,numel(Critical.REA.CCT)) 's \n'];
%     fprintf(strREACCT,Critical.REA.CCT);
%     clear strREACCT

