% =========================================================================
% 教学操作说明：最小梯度点（MGP）搜索或单条轨迹迭代函数
% 使用方法：
%  按以下函数签名调用：function  [deltac_iter,theta_iter,voltage_iter,Normp,no_MGP,flag_normmin]=Fun_Cal_MGP_singletraj_SPM(deltac_start,theta_start,voltage_start,Tunit,n_itermax,norm_Tol,Yfull,preset)
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
function  [deltac_iter,theta_iter,voltage_iter,Normp,no_MGP,flag_normmin]=Fun_Cal_MGP_singletraj_SPM(deltac_start,theta_start,voltage_start,Tunit,n_itermax,norm_Tol,Yfull,preset)
    %% Settings
    G_post=real(Yfull);
    B_post=imag(Yfull);
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    ngen=preset.ngen;
    nbus=preset.nbus;
    %% Initialization
    flag_normmin=0;
    Pe=zeros(n_itermax,ngen);
    Pcoi=zeros(n_itermax,1);
    deltac_iter=zeros(n_itermax,ngen);
    theta_iter=zeros(n_itermax,nbus-ngen);
    voltage_iter=zeros(n_itermax,nbus-ngen);
    Normp=zeros(n_itermax,1);
    norm_min=0;
    no_MGP=0;
    Basevalue = evalin('base','Basevalue');

    %% iteration process
    M = diag([ones(2,1); ones(12,1)*1e-15]);
    options = odeset('Mass',M,'RelTol',1e-10,'AbsTol',[1e-8*ones(1,2),1e-12*ones(1,12)]);
    [theta_start_n,voltage_start_n,flag_iter,n_iter,err] = Fun_AEiteration_SPM(theta_start,voltage_start,deltac_start,preset,Basevalue,"postfault",1e4,1e-10);
    [t_postfault, x_postfault_all] = ode15s(@fred,0:Tunit:Tunit*(n_itermax-1),[deltac_start(2:3); theta_start_n; voltage_start_n],options);
    deltac_temp = x_postfault_all(:,1:2);
    delta1c= -m(2:ngen)'*deltac_temp'/m(1);
    deltac_iter = [delta1c' deltac_temp];
    theta_iter = x_postfault_all(:,3:8);
    voltage_iter = x_postfault_all(:,9:14);



    %% 
    for tm=1:n_itermax
        % Pe calculation
        for  i=1:ngen
            for j=1:ngen
                ddelta=deltac_iter(tm,i)-deltac_iter(tm,j);
                Pe(tm,i)=Pe(tm,i)+E(i)*E(j)*(G_post(i,j)*cos(ddelta)+B_post(i,j)*sin(ddelta));
            end
            for l=1:(nbus-ngen)
                ddelta=deltac_iter(tm,i)-theta_iter(tm,l);
                Pe(tm,i)=Pe(tm,i)+E(i)*voltage_iter(tm,l)*B_post(i,l+ngen)*sin(ddelta)+E(i)*voltage_iter(tm,l)*G_post(i,l+ngen)*cos(ddelta);
            end
        end
        Pcoi(tm)=sum(Pm)-sum(Pe(tm,:));

        
        Normp(tm)=norm((Pm'-Pe(tm,:)-Pcoi(tm,1)/sum(m)*m'));%./d'
        if(tm==1)
            norm_min=Normp(1);
        else
            if(Normp(tm)<norm_min)
                norm_min=Normp(tm);
            end
            if(Normp(tm)-Normp(tm-1)>norm_Tol&&Normp(tm-1)==norm_min&&tm~=2&&norm_min<1e-1)
                flag_normmin=1;
                no_MGP=tm-1;
                fprintf('MGP found since the norm is rather small and a minimum found');
                break;
            end
        end
    end
end
%%
function dfdt = fred(t,x)
    dfdt = f_reducedstate_SPM(x);
end
