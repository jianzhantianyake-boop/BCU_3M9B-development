% =========================================================================
% 教学操作说明：故障轨迹退出点与对应状态的计算函数
% 使用方法：
%  按以下函数签名调用：function [fault_delta_gen,fault_omega,fault_omegac,delta_net_faultclear,voltage_net_faultclear,Exittm,Dotproduct]=Fun_Cal_Exitpoint_SPM(Tlength,Tunit,fault,postfault,preset,delta0,omega0,delta_net0,voltage_net0,Basevalue)
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
function [fault_delta_gen,fault_omega,fault_omegac,delta_net_faultclear,voltage_net_faultclear,Exittm,Dotproduct]=Fun_Cal_Exitpoint_SPM(Tlength,Tunit,fault,postfault,preset,delta0,omega0,delta_net0,voltage_net0,Basevalue)
  %% initialize 
    flag_exit=0;
    cycle=fix(Tlength/Tunit);
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    mT=sum(m,1);
    ngen=preset.ngen;
    nbus=preset.nbus;
    flag_wrap=zeros(ngen,1);
    Pe_post=zeros(cycle,ngen);
    Dotproduct=zeros(cycle,1);
    G_post=real(postfault.Yfull_mod);
    B_post=imag(postfault.Yfull_mod);
    Exittm=0;

%% fault traj
    M = diag([ones(3,1); ones(12,1)*1e-10; ones(3,1)]);
    options = odeset('Mass',M,'RelTol',1e-10,'AbsTol',[1e-8*ones(1,3),1e-12*ones(1,12),1e-8*ones(1,3)]);
    system="fault1";
    delta_net0(fault.faultbus-ngen)=[];
    voltage_net0(fault.faultbus-ngen)=[];
    delta_net0(6)=0;
    voltage_net0(6)=0;
    [delta_net_s,V_net_s,flag_iter,n_iter,err] = Fun_AEiteration_SPM(delta_net0,voltage_net0,delta0,preset,Basevalue,system,1e5,1e-10);
    [t_fault, x_fault_all] = ode15s(@(t,x)f_timedomain_DAE(t, x, system),0:Tunit:Tlength,[delta0; delta_net_s; V_net_s; omega0*ones(ngen,1)],options);
    fault_delta_net=x_fault_all(2:end,4:8);
    fault_delta_net=[fault_delta_net(:,1:(fault.faultbus-ngen-1)) zeros(size(fault_delta_net,1),1) fault_delta_net(:,(fault.faultbus-ngen):end)];
    fault_voltage_net=x_fault_all(2:end,10:14);
    fault_voltage_net=[fault_voltage_net(:,1:(fault.faultbus-ngen-1)) zeros(size(fault_voltage_net,1),1) fault_voltage_net(:,(fault.faultbus-ngen):end)];
    fault_delta_gen = x_fault_all(2:end,1:3);
    fault_omega  = x_fault_all(2:end,16:18);
    fault_omegacoi  = fault_omega *preset.m./sum(preset.m);
    fault_omegac  = fault_omega  - fault_omegacoi *ones(1,3);
    t_fault=t_fault(2:end);
% fault-clear
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
       Results_fsolve=fsolve(@(x)Fun_AEfslove_SPM(x,fault_delta_gen(i,:)',preset,system),temp_ini,optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off','TolX',1e-5));
       delta_net_faultclear(i,:) = Results_fsolve(1:6)';
       voltage_net_faultclear(i,:) = Results_fsolve(7:12)';
       temp_ini=Results_fsolve;
    end   

%% exit point
   for tm=1:cycle
       for  i=1:ngen
            for j=1:ngen
                ddelta=fault_delta_gen(tm,i)-fault_delta_gen(tm,j);
                Pe_post(tm,i)=Pe_post(tm,i)+E(i)*E(j)*(G_post(i,j)*cos(ddelta)+B_post(i,j)*sin(ddelta));
            end
            for l=1:(nbus-ngen)
                ddelta=fault_delta_gen(tm,i)-delta_net_faultclear(tm,l);
                Pe_post(tm,i)=Pe_post(tm,i)+E(i)*voltage_net_faultclear(tm,l)*B_post(i,l+ngen)*sin(ddelta)+E(i)*voltage_net_faultclear(tm,l)*G_post(i,l+ngen)*cos(ddelta);
            end
       end
       Dotproduct(tm)=(Pm'-Pe_post(tm,:))*fault_omegac(tm,:)';
       if tm>2
       if(flag_exit==0 && Dotproduct(tm-1)<0 && Dotproduct(tm)>0)
            flag_exit=1;
            Exittm=tm-1;
       end
       end
   end
end
function dfdt = f_timedomain_DAE(t,x,system)
    dfdt = F_3M9B_SP_DAE(x,system);
end
