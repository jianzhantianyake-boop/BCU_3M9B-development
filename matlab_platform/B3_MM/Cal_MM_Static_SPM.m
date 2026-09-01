% =========================================================================
% 教学操作说明：结构保持模型（SPM）的公共初始化脚本。
% 使用方法：
%   在项目根目录先运行 setup_bcu_paths()，随后由 run_bcu_beginner.m 的
%   "spm_cct"、"spm_numerical" 或 "spm_region" 模式调用，或在已清空的 base workspace 中执行：
%   run(fullfile('B3_MM', 'Cal_MM_Static_SPM.m'))。
% 参数：
%   本文件不接收函数形参。请在下方的设置段配置 preset.PathEnergyCal、
%   preset.EquCal、发电机向量（m、d、Pmpu、xd1、Epu）、ZIP 负荷比例
%   （PloadZIP、QloadZIP）、Case、faultline 和 faultposition。向量长度与
%   发电机顺序必须和所选 case 一致，ZIP 分量的顺序是 Z/I/P。
% 返回 / 工作区结果：
%   在 base workspace 生成 Basevalue、preset、fault、pfdata、netdata、EMF、
%   prefault、fault、postfault 等结构体；其中包含网络节点相角、电压及 DAE
%   后续计算所需的结构保持网络分块数据。
% 步骤：
%   1. 设置能量/平衡点求解策略及机器、ZIP 负荷参数。2. 运行 MATPOWER 潮流
%   并计算内部电势。3. 建立预故障、故障、故障后网络与初始平衡点。
% 单位：
%   角度 rad，角速度基值 rad/s，时间 s；功率、电压、内部电势和导纳通常为 pu。
% 前置条件：
%   需要可用的 MATPOWER、Optimization Toolbox（fsolve）和项目路径；本脚本
%   会 clear 并关闭图窗，后续 SPM CCT/轨迹脚本必须在同一 base workspace 继续运行。
% 研究与验证边界：
%   本次只更新教学注释，不改变 DAE、ZIP 负荷、参数、默认 case 或原始求解逻辑。
%   结果是否收敛及其物理意义必须通过实际 MATLAB 运行和参数记录确认。
% =========================================================================
clear
close all
projectPaths = setup_bcu_paths();
%% 1) 能量计算开关：结构保持模型保留网络节点的能量相关状态。
    preset.PathEnergyCal=20; % 0: Ray approximation  n: N-segment trap approximation (-1)--neglect this part
%% Equilibrium calculation method
    preset.EquCal=2;    % 1--Newton method 2--fsolve method
%% fsolve preset
    options = optimset('TolFun',1e-50,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter','TolX',1e-9);
    options.StepTolerance = 1e-50;
%% 2) 发电机与 ZIP 负荷参数：这里额外保存 Y/I/S 三类负荷分量。
% 9 bus sys===========================================
    Basevalue.omegab=2*pi*60;
    preset.m=[0.1254;0.034;0.016];  %[10;10;10]/omegab;
    preset.d=preset.m.*[0.1;0.1;0.1];
    preset.PloadZIP = [1 0 0]; % Z I P
    preset.QloadZIP = [1 0 0]; % Z I P
% 9 bus sys===========================================
%%%%%%% 思考一下怎么把这几项融到powerflow计算中 %%%%%%%%%
    preset.Pmpu=[0.8980;1.3432;0.9419];
    preset.xd1=[0.0608;0.1198;0.1813];
    preset.Epu=[1.1083;1.1071;1.0606];

    ngen=size(preset.m,1);
    % round(x,3) 表示保留 3 位小数，兼容原先的保留位数意图。
    DHri=round(preset.d./preset.m,3);
    flag_uniform=1;
    for i=2:ngen
        if(preset.d(i)/preset.m(i)~=DHri(1))
            flag_uniform=0;
        end
    end
    preset.flag_uniform=flag_uniform;
    clear flag_uniform DHri
%% 3) MATPOWER 潮流与内部电势：结果将被转换为结构保持模型的初值。
    path_matdata=projectPaths.matpowerData;
    addpath(genpath(path_matdata));
    Case=case9_v2;
    Basevalue.Sbase=Case.baseMVA;
%% run matpower--Matpower powerflow ignore the damping effect, and the results are used to calculate equivalent load
    path_matpower=projectPaths.matpowerRoot;
    addpath(genpath(path_matpower));
    pfdata=Fun_ResultBack(Case);
    if(ngen~=pfdata.bus.numgen)
        error('Powerflow data and generators data not match');
    end
    preset.flagxd=0;    % 0--consider xd' in network already
    % internal EMF calculation
    EMF=Fun_Cal_GenEMF(preset.flagxd,pfdata,preset.xd1);    
    clear path_matpower path_matdata
%% 4) 故障设置：故障母线删除版和大导纳版分别用于不同 DAE 过程。
    preset.faultline=[9;6];  % [Frombus Tobus]
    preset.faultposition=0;
    fault.faultline=preset.faultline;
    fault.faultbus=preset.faultline(preset.faultposition+1,1);
%% 5) 网络构造：同时保存 Yfull、Yfull_forR、Transform 和结构保持分块矩阵。
    %% Add Reactance Xd' into Network
    if(preset.flagxd==0)
        pfdata.branch.RXB_xd=pfdata.branch.RXB;
        fprintf('WARNING: xd already included in casefile!\n');
    else    % add xd' into admittance matrix
        pfdata.branch.RXB_xd=pfdata.branch.RXB;
        for i=1:pfdata.bus.numgen
        no_gen=pfdata.gen.no(i,1);
            for j=1:size(pfdata.branch.RXB_xd,1)
                if(pfdata.branch.RXB_xd(j,1)==no_gen||pfdata.branch.RXB_xd(j,2)==no_gen)
                    pfdata.branch.RXB_xd(j,4)=pfdata.branch.RXB_xd(j,4)+preset.xd1(i);
                end
            end
        end
    end
    %% Transfer RXB into Structure-preserved admittance matrix
    netdata.Y_org=Fun_RXB2Yfull(pfdata.branch.RXB_xd,pfdata);   % admittance without load
    %% Add Passive load into Network
    Yload=zeros(pfdata.bus.numload,5);
    Iload=zeros(pfdata.bus.numload,3);
    Sload=zeros(pfdata.bus.numload,3);
    Yload(:,1)=pfdata.load.no;
    Iload(:,1)=pfdata.load.no;
    Sload(:,1)=pfdata.load.no;
    Y=netdata.Y_org;
    Y_forR=netdata.Y_org;
    for i=1:pfdata.bus.numload
        PLpu=pfdata.load.PQ(i,1)/Basevalue.Sbase;
        QLpu=pfdata.load.PQ(i,2)/Basevalue.Sbase;
        VLpu=pfdata.load.voltage(i,1);
        Yload(i,2)= PLpu*preset.PloadZIP(1)/(VLpu^2);
        Yload(i,3)=-1*QLpu*preset.QloadZIP(1)/(VLpu^2);
        Yload(i,4)= PLpu/(VLpu^2);
        Yload(i,5)=-1*QLpu/(VLpu^2);
        Iload(i,2)= PLpu*preset.PloadZIP(2)/(VLpu);
        Iload(i,3)=-1*QLpu*preset.QloadZIP(2)/(VLpu);
        Sload(i,2)= PLpu*preset.PloadZIP(3);
        Sload(i,3)= QLpu*preset.QloadZIP(3);
    end
    clear PLpu QLpu VLpu
    for i=1:pfdata.bus.numload	% search for load bus
        Loadbus=Yload(i,1);
        for j=1:pfdata.bus.num  %   scan all bus
            if(j==Loadbus)
                Y(j,j)=Y(j,j)+complex(Yload(i,2),Yload(i,3));
                Y_forR(j,j)=Y_forR(j,j)+complex(Yload(i,4),Yload(i,5));
            end
        end
    end
    prefault.Yfull=Y;
    prefault.Yfull_forR=Y_forR;
    preset.Yload = Yload;
    preset.Iload = Iload;
    preset.Sload = Sload;

    preset.genno = pfdata.gen.no;
    
    clear Y Loadbus Y_forR
    % Reduced Network Admittance of Prefault
    [prefault.Yred,prefault.Ynn,prefault.Ynr,prefault.Yrn,prefault.Yrr]=Fun_Yfull2Yred(prefault.Yfull_forR,pfdata,0);
    [prefault.Yfull_mod,prefault.Transform] = Fun_Yfull2Yfull(prefault.Yfull,pfdata,0);
    %% Structure Preserved Admittance of Fault
        Yfull_pre=prefault.Yfull;
        fault.Yfull=Yfull_pre;
        fault.Yfull_mod2 = fault.Yfull;
        fault.Yfull_mod2(fault.faultbus,fault.faultbus)=fault.Yfull_mod2(fault.faultbus,fault.faultbus)+1e6;
        fault.Yfull(:,fault.faultbus)=[];
        fault.Yfull(fault.faultbus,:)=[];

        fault.Yfull_forR=prefault.Yfull_forR;
        fault.Yfull_forR(:,fault.faultbus)=[];
        fault.Yfull_forR(fault.faultbus,:)=[];
        [fault.Yred,fault.Ynn,fault.Ynr,fault.Yrn,fault.Yrr]=Fun_Yfull2Yred(fault.Yfull_forR,pfdata,[1,fault.faultbus]);
        [fault.Yfull_mod,fault.Transform] = Fun_Yfull2Yfull(fault.Yfull_forR,pfdata,[1,fault.faultbus]); % during fault pure impedance %Fun_Yfull2Yfull(fault.Yfull,pfdata,[1,fault.faultbus]);
        [fault.Yfull_mod2,fault.Transform2] = Fun_Yfull2Yfull(fault.Yfull_mod2,pfdata,0);
    %% Structure Preserved Admittance of Postfault
        pfdata.branch.RXB_postfault=pfdata.branch.RXB_xd;
        for i=1:size(pfdata.branch.RXB_postfault,1)
            no1=pfdata.branch.RXB_postfault(i,1);
            no2=pfdata.branch.RXB_postfault(i,2);
            if((fault.faultline(1)==no1&&fault.faultline(2)==no2)||(fault.faultline(1)==no2&&fault.faultline(2)==no1))
                pfdata.branch.RXB_postfault(i,:)=[];
                break;
            end
        end
        clear no1 no2
        postfault.Yorg=Fun_RXB2Yfull(pfdata.branch.RXB_postfault,pfdata);
        Y_forR=postfault.Yorg;
        Y=postfault.Yorg;
        for i=1:pfdata.bus.numload	% search for load bus
            Loadbus=Yload(i,1);
            for j=1:pfdata.bus.num  %   scan all bus
                if(j==Loadbus)
                    Y(j,j)=Y(j,j)+complex(Yload(i,2),Yload(i,3));
                    Y_forR(j,j)=Y_forR(j,j)+complex(Yload(i,4),Yload(i,5));
                end
            end
        end
        postfault.Yfull=Y;
        postfault.Yfull_forR=Y_forR;
        clear Y Loadbus Y_forR
        [postfault.Yfull_mod,postfault.Transform] = Fun_Yfull2Yfull(postfault.Yfull,pfdata,0);
        [postfault.Yred,postfault.Ynn,postfault.Ynr,postfault.Yrn,postfault.Yrr]=Fun_Yfull2Yred(postfault.Yfull_forR,pfdata,0);
        nbus=size(postfault.Yfull,1);
        preset.nbus = nbus;
        preset.ngen = ngen;
%% 6) 结构保持 SEP：未知量包括发电机角度、共同速度、网络角度和电压。
    %% Adopt newton method
    if(preset.EquCal==1)
%         delta0=zeros(ngen,1);
%         [prefault.SEP_delta,prefault.SEP_omegapu,flag_iter,n_iter]=Fun_SEPiteration(prefault.Yfull,preset.Pmpu,preset.Epu,preset.m,preset.d,delta0,0,Basevalue.omegab,1e4,1e-8);
%         if(flag_iter~=1)
%             error(' prefault SEP iteration cannot converage!\n');
%         end
%         clear flag_iter n_iter delta0
%         [postfault.SEP_delta,postfault.SEP_omegapu,flag_iter,n_iter]=Fun_SEPiteration(postfault.Yfull,preset.Pmpu,preset.Epu,preset.m,preset.d,prefault.SEP_delta,prefault.SEP_omegapu-1,Basevalue.omegab,1e4,1e-8);
%         if(flag_iter~=1)
%             error(' postfault SEP iteration cannot converage!\n');
%         end
%         clear flag_iter n_iter
    %% Adopt fsolve
    elseif(preset.EquCal==2)
        x_init=[zeros(ngen,1); zeros((nbus-ngen),1); ones((nbus-ngen),1)];
        Results_fsolve=fsolve(@(x)Fun_SEPfslove_SPM(x,preset,prefault,Basevalue),x_init,options);
        prefault.SEP_omegapu=Results_fsolve(ngen)/Basevalue.omegab+1;
        delta_tmp=[Results_fsolve(1:ngen-1);0];
        deltacoi=delta_tmp'*preset.m/sum(preset.m);
        prefault.SEP_delta=delta_tmp-deltacoi;
        prefault.net_delta=Results_fsolve((ngen+1):nbus)-deltacoi;
        prefault.net_voltage=Results_fsolve((nbus+1):(2*nbus-ngen));
        clear delta_tmp Results_fsolve
        x_init(ngen)=(prefault.SEP_omegapu-1)*Basevalue.omegab;
        for i=1:ngen-1
            x_init(i)=prefault.SEP_delta(i)-prefault.SEP_delta(ngen);
        end
        for i=1:(nbus-ngen)
            x_init(i+ngen)=prefault.net_delta(i)-prefault.SEP_delta(ngen);
        end
        for i=1:(nbus-ngen)
            x_init(i+nbus)=prefault.net_voltage(i);
        end
        Results_fsolve=fsolve(@(x)Fun_SEPfslove_SPM(x,preset,postfault,Basevalue),x_init,options);
        postfault.SEP_omegapu=Results_fsolve(ngen)/Basevalue.omegab+1;
        delta_tmp=[Results_fsolve(1:ngen-1);0];
        deltacoi=delta_tmp'*preset.m/sum(preset.m);
        postfault.SEP_delta=delta_tmp-deltacoi;
        postfault.net_delta=Results_fsolve((ngen+1):nbus)-deltacoi;
        postfault.net_voltage=Results_fsolve((nbus+1):(2*nbus-ngen));
        clear delta_tmp Results_fsolve
    end
%     [prefault.SEP_Perr,flag_SEPerr]=Fun_SEPcheck(prefault,preset,prefault.SEP_delta,(prefault.SEP_omegapu-1)*Basevalue.omegab);
%     if(flag_SEPerr==1)
%         error('SEP calculation error: the SEP solution cannot equilibrium condition! \n')
%     end
%     [postfault.SEP_Perr,flag_SEPerr]=Fun_SEPcheck(postfault,preset,postfault.SEP_delta,(postfault.SEP_omegapu-1)*Basevalue.omegab);
%     if(flag_SEPerr==1)
%         error('SEP calculation error: the SEP solution cannot equilibrium condition! \n')
%     end
    clear flag_SEPerr
    
    clear i j
