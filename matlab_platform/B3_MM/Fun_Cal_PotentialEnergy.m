% =========================================================================
% 教学操作说明：势能分量的计算函数
% 使用方法：
%  按以下函数签名调用：function [Ep1,Ep2,Ep3]=Fun_Cal_PotentialEnergy(preset,postfault,thetac_start,thetac_end)
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
%% 1) 参考项：先计算机械输入与电功率项对应的势能差。
%% Input: thetac_end(n*1)--current point thetac_start--used for path estimation
%% Output: Ep1--equivalent Pm Ep2--equivalent Bij*sinδ item Ep2--path-relevant (lossy network)
function [Ep1,Ep2,Ep3]=Fun_Cal_PotentialEnergy(preset,postfault,thetac_start,thetac_end)
    thetac_SEP=postfault.SEP_delta;
    Y_red=postfault.Yred;
    ngen=size(Y_red,1);
    G_post=real(Y_red);
    B_post=imag(Y_red);
    Pm=preset.Pmpu;
    E=preset.Epu;
    i=preset.m;
    d=preset.d;
    Ep1=0;
    Ep2=0;
    Ep3=0;
    Ep3_rad=0;
    for i=1:ngen
        Ep1=Ep1+(-1)*(Pm(i)-E(i)^2*G_post(i,i))*(thetac_end(i,1)-thetac_SEP(i,1));
    end
    for i=1:ngen-1
        for j=i+1:ngen
            Ep2=Ep2+(-1)*E(i)*E(j)*B_post(i,j)*(cos(thetac_end(i,1)-thetac_end(j,1))-cos(thetac_SEP(i,1)-thetac_SEP(j,1)));
        end
    end
% 2) 路径相关势能：根据当前路径的功角差进行积分。
% radial one-as ref
        for i=1:ngen-1 % 耗散势能
            for j=i+1:ngen
                dtheta_i=thetac_end(i) -thetac_start(i);
                dtheta_j=thetac_end(j) -thetac_start(j);
                dtheta_ij=dtheta_i-dtheta_j;
                if(abs(dtheta_ij)>1e-6) % 避免数值病态问题
                    adb=(dtheta_i+dtheta_j)/dtheta_ij;                
                else
                    adb=dtheta_i+dtheta_j;
                end
                Ep3_rad= Ep3_rad+E(i)*E(j)*G_post(i,j)*adb*(sin(thetac_end(i)-thetac_end(j))-sin(thetac_start(i)-thetac_start(j)));
            end
        end



    % 3) PathEnergyCal=0 使用 Ray 近似；正整数使用中间点梯形积分。
    if(preset.PathEnergyCal==0)
        for i=1:ngen-1 % 耗散势能
            for j=i+1:ngen
                dtheta_i=thetac_end(i) -thetac_start(i);
                dtheta_j=thetac_end(j) -thetac_start(j);
                dtheta_ij=dtheta_i-dtheta_j;
                if(abs(dtheta_ij)>1e-7) % 避免数值病态问题
                    adb=(dtheta_i+dtheta_j)/dtheta_ij;                
                else
                    adb=dtheta_i+dtheta_j;
                end
                Ep3= Ep3+E(i)*E(j)*G_post(i,j)*adb*(sin(thetac_end(i)-thetac_end(j))-sin(thetac_start(i)-thetac_start(j)));
            end
        end
    elseif(preset.PathEnergyCal==-1)
        Ep3=0;
    else
        n_mid=preset.PathEnergyCal-1;   % numbers of inserted mid-point (n_mid=0--trapezoidal from start to end)  
        dtheta=thetac_end-thetac_start;  % dtheta(i)=theta_end(i)-theta_start(i)
        unit_dtheta=dtheta/(n_mid+1);
        for i=1:ngen
            for j=1:ngen
                if(i~=j)
                    for m=1:n_mid+1
                    Ep3=Ep3+E(i)*E(j)*G_post(i,j)*0.5*unit_dtheta(i)...,
                        *(cos(thetac_start(i)+(m-1)*unit_dtheta(i)-thetac_start(j)-(m-1)*unit_dtheta(j))+...,
                        cos(thetac_start(i)+m*unit_dtheta(i)-thetac_start(j)-m*unit_dtheta(j)));
                    end
                end
            end
        end
    end

end
