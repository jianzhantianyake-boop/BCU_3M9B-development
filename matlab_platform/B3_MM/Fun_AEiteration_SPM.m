% =========================================================================
% 教学操作说明：结构保持模型代数方程的残差或迭代求解函数
% 使用方法：
%  按以下函数签名调用：function [delta_net_s,V_net_s,flag_iter,n_iter,err] = Fun_AEiteration_SPM(delta_net0,V_net0,deltac,preset,basevalue,system,n_itermax,Tolerr)
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
% x(2*nbus-ngen) = delta_net(1:nbus-ngen) | V_net(nbus-ngen+1:2*nbus-2*ngen)
function [delta_net_s,V_net_s,flag_iter,n_iter,err] = Fun_AEiteration_SPM(delta_net0,V_net0,deltac,preset,basevalue,system,n_itermax,Tolerr)
%% parameter preprocess
    %system = evalin('base','system');
    prefault = evalin('base','prefault');
    fault = evalin('base','fault');
    postfault = evalin('base','postfault');
    switch system
        case "prefault"
            Yfull = prefault.Yfull_mod;
            Transform = prefault.Transform;
        case "fault1"
            Yfull = fault.Yfull_mod;
            Transform = fault.Transform;
        case "fault2"
            Yfull = fault.Yfull_mod2;
            Transform = fault.Transform2;
        case "postfault"
            Yfull = postfault.Yfull_mod;
            Transform = postfault.Transform;
    end
    G=real(Yfull);    
    B=imag(Yfull);
    ngen=size(preset.genno,1);
    nbus=size(Yfull,1);
    E=preset.Epu;

%% Initialization
    V_net=V_net0;
    delta_net=delta_net0;

%% iteration settings
    n_iter=1;
    flag_iter=0;    % 0--uncompleted 1--success 2--failed(n_iter>n_itermax)
    Step_len=1;
%% iteration
while(flag_iter==0)
    Pnet = zeros(nbus-ngen,1);
    Qnet = zeros(nbus-ngen,1);
    % P calculation of Bus
    for i=1:(nbus-ngen)
        for j=1:ngen
            ddelta=delta_net(i)-deltac(j);
            Pnet(i)=Pnet(i)+V_net(i)*E(j)*B(i+ngen,j)*sin(ddelta)+V_net(i)*E(j)*G(i+ngen,j)*cos(ddelta);
        end
        for l=1:(nbus-ngen)
            ddelta=delta_net(i)-delta_net(l);
            Pnet(i)=Pnet(i)+V_net(i)*V_net(l)*B(i+ngen,l+ngen)*sin(ddelta)+V_net(i)*V_net(l)*G(i+ngen,l+ngen)*cos(ddelta);
        end
        for h=1:size(preset.Sload,1)
           if (preset.Sload(h,1)==Transform(i+ngen))
               if (system == "fault1")||(system == "fault2")
                   if Transform(i+ngen)~=fault.faultbus
                        %Pnet(i)=Pnet(i)+preset.Sload(h,2); % during fault
                        %pure impedance
                   end
               else
                   Pnet(i)=Pnet(i)+preset.Sload(h,2);
               end
           end
        end
    end
    % Q calculation of Bus
    for i=1:(nbus-ngen)
        for j=1:ngen
            ddelta=delta_net(i)-deltac(j);
            Qnet(i)=Qnet(i)-V_net(i)*E(j)*B(i+ngen,j)*cos(ddelta)+V_net(i)*E(j)*G(i+ngen,j)*sin(ddelta);
        end
        for l=1:(nbus-ngen)
            ddelta=delta_net(i)-delta_net(l);
            Qnet(i)=Qnet(i)-V_net(i)*V_net(l)*B(i+ngen,l+ngen)*cos(ddelta)+V_net(i)*V_net(l)*G(i+ngen,l+ngen)*sin(ddelta);
        end
        for h=1:size(preset.Sload,1)
            if (preset.Sload(h,1)==Transform(i+ngen))
               if (system == "fault1")||(system == "fault2")
                   if Transform(i+ngen)~=fault.faultbus
                        %Qnet(i)=Qnet(i)+preset.Sload(h,3); % during fault
                        %pure impedance
                   end
               else
                   Qnet(i)=Qnet(i)+preset.Sload(h,3);
               end
               
            end
        end
    end
    % Jacobi matrix calculation
    J11=zeros((nbus-ngen),(nbus-ngen)); % Pnet/delta
    J12=zeros((nbus-ngen),(nbus-ngen)); % Pnet/V
    J21=zeros((nbus-ngen),(nbus-ngen)); % Qnet/delta
    J22=zeros((nbus-ngen),(nbus-ngen)); % Qnet/V
    for i=1:(nbus-ngen)
        J12(i,i)=2*V_net(i)*G(i+ngen,i+ngen);
        J22(i,i)=-2*V_net(i)*B(i+ngen,i+ngen);
        for j=1:ngen
            ddelta=delta_net(i)-deltac(j);
            J11(i,i)=J11(i,i)+V_net(i)*E(j)*B(i+ngen,j)*cos(ddelta)-V_net(i)*E(j)*G(i+ngen,j)*sin(ddelta);
            J12(i,i)=J12(i,i)+E(j)*B(i+ngen,j)*sin(ddelta)+E(j)*G(i+ngen,j)*cos(ddelta);
            J21(i,i)=J21(i,i)+V_net(i)*E(j)*B(i+ngen,j)*sin(ddelta)+V_net(i)*E(j)*G(i+ngen,j)*cos(ddelta);
            J22(i,i)=J22(i,i)-E(j)*B(i+ngen,j)*cos(ddelta)+E(j)*G(i+ngen,j)*sin(ddelta);
        end
        for j=1:(nbus-ngen)
            if (i~=j)
                ddelta=delta_net(i)-delta_net(j);
                J11(i,i)=J11(i,i)+V_net(i)*V_net(j)*B(i+ngen,j+ngen)*cos(ddelta)-V_net(i)*V_net(j)*G(i+ngen,j+ngen)*sin(ddelta);
                J11(i,j)=-V_net(i)*V_net(j)*B(i+ngen,j+ngen)*cos(ddelta)+V_net(i)*V_net(j)*G(i+ngen,j+ngen)*sin(ddelta);

                J12(i,i)=J12(i,i)+V_net(j)*B(i+ngen,j+ngen)*sin(ddelta)+V_net(j)*G(i+ngen,j+ngen)*cos(ddelta);
                J12(i,j)=V_net(i)*B(i+ngen,j+ngen)*sin(ddelta)+V_net(i)*G(i+ngen,j+ngen)*cos(ddelta);

                J21(i,i)=J21(i,i)+V_net(i)*V_net(j)*B(i+ngen,j+ngen)*sin(ddelta)+V_net(i)*V_net(j)*G(i+ngen,j+ngen)*cos(ddelta);
                J21(i,j)=-V_net(i)*V_net(j)*B(i+ngen,j+ngen)*sin(ddelta)-V_net(i)*V_net(j)*G(i+ngen,j+ngen)*cos(ddelta);

                J22(i,i)=J22(i,i)-V_net(j)*B(i+ngen,j+ngen)*cos(ddelta)+V_net(j)*G(i+ngen,j+ngen)*sin(ddelta);
                J22(i,j)=-V_net(i)*B(i+ngen,j+ngen)*cos(ddelta)+V_net(i)*G(i+ngen,j+ngen)*sin(ddelta);
            end
        end
    end
    J=[J11, J12;J21, J22];
    dx= -1* J \ [Pnet ; Qnet]*Step_len;
    delta_net(1:(nbus-ngen)) = delta_net(1:(nbus-ngen))+dx(1:(nbus-ngen));
    V_net(1:(nbus-ngen)) = V_net(1:(nbus-ngen)) + dx((nbus-ngen)+1:2*(nbus-ngen));
    err = max(abs([Pnet ; Qnet]));

    if(abs(err)<Tolerr)
        flag_iter=1;
    else
        n_iter=n_iter+1;
    end
    if(n_iter>=n_itermax)
        flag_iter=2;
    end

end

if(flag_iter==1)
    delta_net_s=delta_net;
    V_net_s=V_net;
else
    delta_net_s=nan;
    V_net_s=nan;
end


end
