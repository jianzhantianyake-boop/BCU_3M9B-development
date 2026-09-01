% =========================================================================
% 教学操作说明：网络约简模型的公共初始化脚本。
% 使用方法：
%   在项目根目录先运行 setup_bcu_paths()，随后由 run_bcu_beginner.m 的
%   "reduced_cct"、"reduced_numerical" 或 "reduced_region" 模式调用，或在已清空的 base workspace 中执行：
%   run(fullfile('B3_MM', 'Cal_MM_Static.m'))。
% 参数：
%   本文件不接收函数形参。请在下方的设置段依次配置 preset.PathEnergyCal、
%   preset.EquCal、preset.m、preset.d、preset.Pmpu、preset.xd1、preset.Epu、
%   Case、preset.faultline 与 preset.faultposition。五个发电机向量必须与
%   所选 case 中发电机编号顺序和数量一致。
% 返回 / 工作区结果：
%   在 base workspace 生成 Basevalue、preset、fault、pfdata、netdata、EMF、
%   prefault、fault、postfault 等结构体，以及预故障/故障/故障后网络与 SEP。
% 步骤：
%   1. 配置能量计算、平衡点求解器和发电机参数。2. 运行 MATPOWER 潮流并
%   计算内部电势。3. 构造三种网络并求解 SEP，供 CCT、轨迹或区域搜索使用。
% 单位：
%   角度 rad，角速度基值 rad/s，时间 s；功率、电压、内部电势和导纳通常为 pu。
% 前置条件：
%   需要可用的 MATPOWER、Optimization Toolbox（fsolve）和项目路径；本脚本
%   会关闭图窗，且其结果必须保留在 base workspace，不能在函数局部调用后丢失。
% 研究与验证边界：
%   本次仅把注释改为操作说明，不改变任何初始化方程、默认参数、故障定义或求解器。
%   切换到 39 母线前必须完成向量维度与结果的实际 MATLAB 验证。
% =========================================================================
close all
projectPaths = setup_bcu_paths();
BCUCFG = bcu_config();   % 集中配置（参数覆盖来源；缺字段时用下方原始默认，零回归）
%% 1) 能量计算开关：这里只设置后续能量函数的路径积分策略。
% 0 表示 Ray approximation，正整数表示分段梯形积分，-1 表示省略该项。
    preset.PathEnergyCal=bcu_pick(BCUCFG,'PathEnergyCal',0); % 0: Ray approximation  n: N-segment trap approximation (-1)--neglect this part
%% Equilibrium calculation method
    preset.EquCal=bcu_pick(BCUCFG,'EquCal',2);    % 1--Newton method 2--fsolve method
%% fsolve preset
    options = optimset('TolFun',bcu_pick(BCUCFG,'TolFun',1e-50),'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter','TolX',bcu_pick(BCUCFG,'TolX',1e-9));
    options.StepTolerance = 1e-50;
%% 2) 发电机与基值参数：m、d、Pm、xd1、Epu 必须与发电机顺序一致。
% 9 bus sys===========================================
    Basevalue.omegab=2*pi*bcu_pick(BCUCFG,'f_base',60);
    preset.m=bcu_pick(BCUCFG,'m',[0.1254;0.034;0.016]);  %[10;10;10]/omegab; preset.m=[0.1254;0.034;0.016];[0.016;0.034;0.1254]
    preset.d=preset.m.*bcu_pick(BCUCFG,'damping_ratio',[0.1;0.1;0.1]);
% % 39 bus sys=========================================
%     Basevalue.omegab=2*pi*50;
%     preset.m=[42.0;30.3;35.8;28.6;26;34.8;26.4;24.3;34.5;50]*2/(Basevalue.omegab);
%     preset.d=preset.m*0.05;
%     no_GFM=[3;4;10];
%     for i=1:size(no_GFM,1)
%         preset.d(no_GFM(i))=preset.m(no_GFM(i))*0.8;
% %     preset.d(no_GFM(i))=50/Basevalue.omegab;
%     end
% %     preset.d(2)=preset.m(2)*0.7;
% %     preset.d(2)=50/Basevalue.omegab;
%     clear no_GFM
%     
%     DT=sum(preset.d,1);
%     HT=sum(preset.m,1);
%%%%%%% 思考一下怎么把这几项融到powerflow计算中 %%%%%%%%%
% 9 bus sys===========================================
    preset.Pmpu=bcu_pick(BCUCFG,'Pm',[0.8980;1.3432;0.9419]);
    preset.xd1=bcu_pick(BCUCFG,'xd1',[0.0608;0.1198;0.1813]);
    preset.Epu=bcu_pick(BCUCFG,'E',[1.1083;1.1071;1.0606]);

% 39 bus sys=========================================
% order: bus no.30 to bus no.39
%     preset.Pmpu=[2.50;6.77871;6.5;6.32;5.08;6.5;5.6;5.4;8.3;10];
%     preset.xd1=[0.031;0.0697;0.0531;0.0436;0.132;0.05;0.049;0.057;0.057;0.006];
%     preset.Epu=[1.0929;1.1966;1.1491;1.0808;1.3971;1.1910;1.1394;1.0709;1.1368;1.0368];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% 3) MATPOWER 潮流：先求稳态电压/功率，再把结果转换成 BCU 的 pfdata。
    path_matdata=projectPaths.matpowerData;
    addpath(genpath(path_matdata));
    Case=feval(bcu_pick(BCUCFG,'CaseName','case9_v2'));
%     Case=case39_modified;
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
%% 4) 故障设置：faultline 给出线路两端，faultposition 决定故障母线。
%     preset.flagfault=1; % 1--fault occur
    preset.faultline=bcu_pick(BCUCFG,'faultline',[9;6]);  % [Frombus Tobus]
    preset.faultposition=bcu_pick(BCUCFG,'faultposition',0);
    fault.faultline=preset.faultline;
%     fault.flagfault=preset.flagfault; % 1--fault occur
    fault.faultbus=preset.faultline(preset.faultposition+1,1);
%% 5) 网络构造：加入负荷和暂态电抗，分别形成预故障、故障、故障后网络。
% 后续 Yred 的维度和节点顺序都由这一段决定。
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
    Yload=zeros(pfdata.bus.numload,3);
    Yload(:,1)=pfdata.load.no;
    Y=netdata.Y_org;
    for i=1:pfdata.bus.numload
        PLpu=pfdata.load.PQ(i,1)/Basevalue.Sbase;
        QLpu=pfdata.load.PQ(i,2)/Basevalue.Sbase;
        VLpu=pfdata.load.voltage(i,1);
        Yload(i,2)=PLpu/(VLpu^2);
        Yload(i,3)=-1*QLpu/(VLpu^2);
    end
    clear PLpu QLpu VLpu
    for i=1:pfdata.bus.numload	% search for load bus
        Loadbus=Yload(i,1);
        for j=1:pfdata.bus.num  %   scan all bus
            if(j==Loadbus)
                Y(j,j)=Y(j,j)+complex(Yload(i,2),Yload(i,3));
            end
        end
    end
    prefault.Yfull=Y;
    
    clear Y Loadbus
    %% Reduced Network Admittance of Prefault
        [prefault.Yred,prefault.Ynn,prefault.Ynr,prefault.Yrn,prefault.Yrr]=Fun_Yfull2Yred(prefault.Yfull,pfdata,0);
    %% Reduced Network Admittance of Fault
        Yfull_pre=prefault.Yfull;
        fault.Yfull=Yfull_pre;
        fault.Yfull(:,fault.faultbus)=[];
        fault.Yfull(fault.faultbus,:)=[];
        [fault.Yred,fault.Ynn,fault.Ynr,fault.Yrn,fault.Yrr]=Fun_Yfull2Yred(fault.Yfull,pfdata,[1,fault.faultbus]);
    %% Reduced Network Admittance of Postfault
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
        Y=postfault.Yorg;
        for i=1:pfdata.bus.numload	% search for load bus
        Loadbus=Yload(i,1);
        for j=1:pfdata.bus.num  %   scan all bus
            if(j==Loadbus)
                Y(j,j)=Y(j,j)+complex(Yload(i,2),Yload(i,3));
            end
        end
        end
        postfault.Yfull=Y;
        clear Y Loadbus
        [postfault.Yred,postfault.Ynn,postfault.Ynr,postfault.Yrn,postfault.Yrr]=Fun_Yfull2Yred(postfault.Yfull,pfdata,0);
%% 6) SEP 求解：用 Newton 或 fsolve 求预故障/故障后平衡点，并检查残差。
    %% Adopt newton method
    if(preset.EquCal==1)
        delta0=zeros(ngen,1);
        [prefault.SEP_delta,prefault.SEP_omegapu,flag_iter,n_iter]=Fun_SEPiteration(prefault.Yred,preset.Pmpu,preset.Epu,preset.m,preset.d,delta0,0,Basevalue.omegab,1e4,1e-8);
        if(flag_iter~=1)
            error(' prefault SEP iteration cannot converage!\n');
        end
        %clear flag_iter n_iter delta0
        [postfault.SEP_delta,postfault.SEP_omegapu,flag_iter,n_iter]=Fun_SEPiteration(postfault.Yred,preset.Pmpu,preset.Epu,preset.m,preset.d,prefault.SEP_delta,prefault.SEP_omegapu-1,Basevalue.omegab,1e4,1e-8);
        if(flag_iter~=1)
            error(' postfault SEP iteration cannot converage!\n');
        end
        %clear flag_iter n_iter
    %% Adopt fsolve
    elseif(preset.EquCal==2)
        deltaomega_init=zeros(ngen,1);
        Results_fsolve=fsolve(@(delta_omega)Fun_SEPfslove(delta_omega,preset,prefault,Basevalue),deltaomega_init,options);
        prefault.SEP_omegapu=Results_fsolve(ngen)/Basevalue.omegab+1;
        delta_tmp=[Results_fsolve(1:ngen-1);0];
        prefault.SEP_delta=delta_tmp-delta_tmp'*preset.m/sum(preset.m);
        clear delta_tmp Results_fsolve
        deltaomega_init(ngen)=(prefault.SEP_omegapu-1)*Basevalue.omegab;
        for i=1:ngen-1
            deltaomega_init(i)=prefault.SEP_delta(i)-prefault.SEP_delta(ngen);
        end
        Results_fsolve=fsolve(@(delta_omega)Fun_SEPfslove(delta_omega,preset,postfault,Basevalue),deltaomega_init,options);
        postfault.SEP_omegapu=Results_fsolve(ngen)/Basevalue.omegab+1;
        delta_tmp=[Results_fsolve(1:ngen-1);0];
        postfault.SEP_delta=delta_tmp-delta_tmp'*preset.m/sum(preset.m);
        clear delta_tmp Results_fsolve deltaomega_init
    end
    [prefault.SEP_Perr,flag_SEPerr]=Fun_SEPcheck(prefault,preset,prefault.SEP_delta,(prefault.SEP_omegapu-1)*Basevalue.omegab);
    if(flag_SEPerr==1)
        error('SEP calculation error: the SEP solution cannot equilibrium condition! \n')
    end
    [postfault.SEP_Perr,flag_SEPerr]=Fun_SEPcheck(postfault,preset,postfault.SEP_delta,(postfault.SEP_omegapu-1)*Basevalue.omegab);
    if(flag_SEPerr==1)
        error('SEP calculation error: the SEP solution cannot equilibrium condition! \n')
    end
    clear flag_SEPerr
    
    clear i j
