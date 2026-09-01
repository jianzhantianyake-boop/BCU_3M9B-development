% =========================================================================
% 教学操作说明：势能分量的计算函数
% 使用方法：
%  按以下函数签名调用：function [Ep1,Ep2,Ep3,Ep4,Ep5]=Fun_Cal_PotentialEnergy_SPM(preset,postfault,thetac_end,theta_net_end,voltage_net_end)
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
%% 1) 结构保持势能：同时考虑发电机、网络节点和负荷相关项。
%% Input: thetac_end(n*1)--current point thetac_start--used for path estimation
%% Output: Ep1--generation potential EP2--network potential Ep3-load losses
function [Ep1,Ep2,Ep3,Ep4,Ep5]=Fun_Cal_PotentialEnergy_SPM(preset,postfault,thetac_end,theta_net_end,voltage_net_end)
    thetac_SEP=postfault.SEP_delta;
    theta_net_SEP=postfault.net_delta;
    voltage_net_SEP=postfault.net_voltage;
    Yfull = postfault.Yfull_mod;
    G_post=real(Yfull);
    B_post=imag(Yfull);
    ngen=size(preset.genno,1);
    nbus=size(Yfull,1);
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    Ep1=0;
    Ep2=0;
    Ep3=0;
    Ep4=0;
    Ep5=0;


    for i=1:ngen
        Ep1=Ep1+(-1)*(Pm(i)-E(i)^2*G_post(i,i))*(thetac_end(i,1)-thetac_SEP(i,1));
    end
    for i=1:(nbus-ngen)
        for j=1:ngen
            Ep2=Ep2+(-1)*(voltage_net_end(i,1)*E(j)*B_post(i+ngen,j)*cos(theta_net_end(i,1)-thetac_end(j,1))-voltage_net_SEP(i,1)*E(j)*B_post(i+ngen,j)*cos(theta_net_SEP(i,1)-thetac_SEP(j,1)) );
        end
        Ep2=Ep2+(-1)*(voltage_net_end(i,1)^2/2*B_post(i+ngen,i+ngen)-voltage_net_SEP(i,1)^2/2*B_post(i+ngen,i+ngen));
        for h=1:size(preset.Sload,1)
            if (preset.Sload(h,1)==postfault.Transform(i+ngen))
                 Ep5=Ep5+preset.Sload(h,2)*(theta_net_end(i,1)-theta_net_SEP(i,1))+preset.Sload(h,3)*(log(voltage_net_end(i,1))- log(voltage_net_SEP(i,1)));
            end
        end
    end
    for i=1:(nbus-ngen-1)
        for l=i+1:(nbus-ngen)
            Ep2=Ep2+(-1)*(voltage_net_end(i,1)*voltage_net_end(l,1)*B_post(i+ngen,l+ngen)*cos(theta_net_end(i,1)-theta_net_end(l,1))-voltage_net_SEP(i,1)*voltage_net_SEP(l,1)*B_post(i+ngen,l+ngen)*cos(theta_net_SEP(i,1)-theta_net_SEP(l,1)) );
        end
    end

    for i=1:(nbus-ngen)
         Ep3=Ep3+G_post(i+ngen,i+ngen)/3*(theta_net_end(i)-theta_net_SEP(i))*( voltage_net_end(i)^2 + voltage_net_end(i)*voltage_net_SEP(i) + voltage_net_SEP(i)^2 );
    end


    % network losses
    % 2) 路径策略：Ray 近似与分段积分只改变数值近似，不改变状态模型。
    if(preset.PathEnergyCal==0)
        Ep4=0;  
        %Ray approximation: to be completed
    elseif(preset.PathEnergyCal==-1)
        Ep4=0;
    else
        n_mid=preset.PathEnergyCal-1;   % numbers of inserted mid-point (n_mid=0--trapezoidal from start to end)  
        dtheta=thetac_end-thetac_SEP;  % dtheta(i)=theta_end(i)-theta_start(i)
        dtheta_net=theta_net_end-theta_net_SEP;
        dvoltage=voltage_net_end-voltage_net_SEP;
        unit_dtheta=dtheta/(n_mid+1);
        unit_dtheta_net=dtheta_net/(n_mid+1);
        unit_dvoltage=dvoltage/(n_mid+1);
        % P lossy of gen
        for i=1:ngen
            for j=1:ngen
                if(i~=j)
                    for m=1:n_mid+1

                        Ep4=Ep4+E(i)*E(j)*G_post(i,j)*0.5*unit_dtheta(i)...,
                            *(cos(thetac_SEP(i)+(m-1)*unit_dtheta(i)-thetac_SEP(j)-(m-1)*unit_dtheta(j))...,
                            +cos(thetac_SEP(i)+m*unit_dtheta(i)-thetac_SEP(j)-m*unit_dtheta(j)));
                    end
                end
            end
            for l=1:(nbus-ngen)
                for m=1:n_mid+1
                    Ep4=Ep4+E(i)*G_post(i,l+ngen)*0.5*unit_dtheta(i)...,
                        *((voltage_net_SEP(l)+(m-1)*unit_dvoltage(l))*cos(thetac_SEP(i)+(m-1)*unit_dtheta(i)-theta_net_SEP(l)-(m-1)*unit_dtheta_net(l))...,
                            +(voltage_net_SEP(l)+m*unit_dvoltage(l))*cos(thetac_SEP(i)+m*unit_dtheta(i)-theta_net_SEP(l)-m*unit_dtheta_net(l)));
                end
            end 
        end
        % P lossy of Bus
        for i=1:(nbus-ngen)
            for j=1:ngen
                for m=1:n_mid+1
                    Ep4=Ep4+E(j)*G_post(i+ngen,j)*0.5*unit_dtheta_net(i)...,
                        *((voltage_net_SEP(i)+(m-1)*unit_dvoltage(i))*cos(theta_net_SEP(i)+(m-1)*unit_dtheta_net(i)-thetac_SEP(j)-(m-1)*unit_dtheta(j))...,
                        +(voltage_net_SEP(i)+m*unit_dvoltage(i))*cos(theta_net_SEP(i)+m*unit_dtheta_net(i)-thetac_SEP(j)-m*unit_dtheta(j)));
                end
            end
            for l=1:(nbus-ngen)
                if(i~=l)
                    for m=1:n_mid+1
                        Ep4=Ep4+G_post(i+ngen,l+ngen)*0.5*unit_dtheta_net(i)...,
                        *((voltage_net_SEP(i)+(m-1)*unit_dvoltage(i))*(voltage_net_SEP(l)+(m-1)*unit_dvoltage(l))*cos(theta_net_SEP(i)+(m-1)*unit_dtheta_net(i)-theta_net_SEP(l)-(m-1)*unit_dtheta_net(l))...,
                        +(voltage_net_SEP(i)+(m)*unit_dvoltage(i))*(voltage_net_SEP(l)+(m)*unit_dvoltage(l))*cos(theta_net_SEP(i)+(m)*unit_dtheta_net(i)-theta_net_SEP(l)-(m)*unit_dtheta_net(l)));
                    end
                end
            end
        end
        % Q/V lossy of Bus
        for i=1:(nbus-ngen)
            for j=1:ngen
                for m=1:n_mid+1
                    Ep4 = Ep4+E(j)*G_post(i+ngen,j)*0.5*unit_dvoltage(i)...,
                        *( sin(theta_net_SEP(i)+(m-1)*unit_dtheta_net(i)-thetac_SEP(j)-(m-1)*unit_dtheta(j))...,
                        + sin(theta_net_SEP(i)+(m)*unit_dtheta_net(i)-thetac_SEP(j)-(m)*unit_dtheta(j)) );
                end
            end
            for l=1:(nbus-ngen)
                if(i~=l)
                    for m=1:n_mid+1
                        Ep4 = Ep4+G_post(i+ngen,l+ngen)*0.5*unit_dvoltage(i)...,
                            *((voltage_net_SEP(l)+(m-1)*unit_dvoltage(l))*sin(theta_net_SEP(i)+(m-1)*unit_dtheta_net(i)-theta_net_SEP(l)-(m-1)*unit_dtheta_net(l)) ...,
                            +(voltage_net_SEP(l)+m*unit_dvoltage(l))*sin(theta_net_SEP(i)+m*unit_dtheta_net(i)-theta_net_SEP(l)-m*unit_dtheta_net(l))  );
                    end
                end
            end
        end  

    end

end
