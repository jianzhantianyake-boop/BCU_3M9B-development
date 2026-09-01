% =========================================================================
% 教学操作说明：稳定平衡点的迭代计算函数
% 使用方法：
%  按以下函数签名调用：function [SEP,omegacoi,flag_iter,n_iter]=Fun_SEPiteration_SPM(Y_full,pfdata,Pm,E,m,d,delta0,omega0,omegab,n_itermax,Tolerr)
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
function [SEP,omegacoi,flag_iter,n_iter]=Fun_SEPiteration_SPM(Y_full,pfdata,Pm,E,m,d,delta0,omega0,omegab,n_itermax,Tolerr)
%% parameters preprocess  
    G=real(Y_full);    
    B=imag(Y_full);
    ngen=size(Y_red,1);
    HT=sum(m,1);
    DT=sum(d,1);
%% initialize
    deltan=delta0(ngen,1);
    thetatt=zeros(ngen,1);
    Step_len=5e-2;
    n_cntmin=0;
    Perr=inf;
% calculate Pcoi
    Pett=0;
    for i=1:ngen
        Pett=Pett+G(i,i)*E(i)^2;
        for j=1:ngen
            if(i~=j)
                Pett=Pett+E(i)*E(j)*G(i,j)*cos(delta0(i)-delta0(j))...,
                    +E(i)*E(j)*B(i,j)*sin(delta0(i)-delta0(j));
            end
        end
    end
    Pcoi=sum(Pm)-Pett;
% Initial Pei
    Pe_ref0=Pm-Pcoi/DT*d;
    delta=delta0(1:ngen-1,1);   % (n-1) delta
    omegacoi=omega0;
% iteration settings
    n_iter=1;
    flag_iter=0;    % 0--uncompleted 1--success 2--failed(n_iter>n_itermax)
%% Powerflow calculation--Reduced matrix
while(flag_iter==0)
        Pe_ref=Pe_ref0(1:ngen-1,1);
    % fi(delta_k(1),...,delta_k(n-1)) calculation
        f0=zeros(ngen-1,1);
        for i=1:ngen-1
            f0(i)=E(i)^2*G(i,i)+E(i)*E(ngen)*G(i,ngen)*cos(delta(i)-deltan)...,
                 +E(i)*E(ngen)*B(i,ngen)*sin(delta(i)-deltan);
            for j=1:ngen-1
                if(i~=j)
                    f0(i)=f0(i)+E(i)*E(j)*G(i,j)*cos(delta(i)-delta(j))...,
                         +E(i)*E(j)*B(i,j)*sin(delta(i)-delta(j));
                end
            end
        end
    % Jacobi matrix calculation
        J=zeros(ngen-1,ngen-1);
        for i=1:ngen-1
            J(i,i)=E(i)*E(ngen)*(-G(i,ngen)*sin(delta(i)-deltan)+B(i,ngen)*cos(delta(i)-deltan));
            for j=1:ngen-1
                if (i~=j)
                    J(i,i)=J(i,i)+E(i)*E(j)*(-G(i,j)*sin(delta(i)-delta(j))+B(i,j)*cos(delta(i)-delta(j)));
                    J(i,j)=E(i)*E(j)*(G(i,j)*sin(delta(i)-delta(j))-B(i,j)*cos(delta(i)-delta(j)));
                end
            end
        end
        delta_his=delta;
        deltan_his=deltan;
        ddelta=inv(J)*(Pe_ref-f0)*Step_len;
        delta=delta+ddelta;
        thetatt(1:ngen-1,1)=delta;
        thetatt(ngen,1)=deltan;
        theta_coi=m'*thetatt/HT;
    %     if(theta_coi>2*pi)
    %         thetatt(:,1)=thetatt(:,1)-2*pi*ones(ngen,1);
    %         theta_coi=theta_coi-2*pi;
    %     end
        deltatt=thetatt-theta_coi*ones(ngen,1);
        deltan=deltatt(ngen,1);
        delta=deltatt(1:ngen-1,1);
    % calculate Pe(n)
        Pen=E(ngen)^2*G(ngen,ngen);
        for j=1:ngen-1
            Pen=Pen+E(ngen)*E(j)*G(ngen,j)*cos(deltan-delta(j))...,
                +E(ngen)*E(j)*B(ngen,j)*sin(deltan-delta(j));
        end
    % check for error
        Perr_his=Perr;
        Perr=Pen-Pe_ref0(ngen); % the additional power inject n
%         if(abs(Perr)>abs(Perr_his))
%             Step_len=Step_len/1.1;
%             if(Step_len<=1e-5)
%                 Step_len=1e-5;
%                 n_cntmin=n_cntmin+1;
%             end
%             flag_repeat=1;
%             delta=delta_his;
%             deltan=deltan_his;
%             Perr=Perr_his;
%         else
%             flag_repeat=0;
%         end
%         if(flag_repeat==1&&Step_len==1e-5&&n_cntmin>2)
%             error('Iteration failed: Perr cannot further reduced with the minimized step_len\n');
%         end
    flag_repeat=0;
    % update P
    if(flag_repeat~=1)
    % recalculate Pcoi
        Pett=0;
        for i=1:ngen
            Pett=Pett+G(i,i)*E(i)^2;
            for j=1:ngen
                if(i~=j)
                    Pett=Pett+E(i)*E(j)*G(i,j)*cos(thetatt(i)-thetatt(j))...,
                        +E(i)*E(j)*B(i,j)*sin(thetatt(i)-thetatt(j));
                end                   
            end
        end
        Pcoi=sum(Pm)-Pett;
        % whether error is sufficiently small 
        if(abs(Perr)<Tolerr)
            flag_iter=1;
        else
            n_iter=n_iter+1;        
            Pe_ref0=Pm-(Pcoi-Perr)/sum(d)*d;
        end
        if(n_iter>=n_itermax)
            flag_iter=2;
        end
    end
end

    if(flag_iter==1)
        SEP=deltatt;
        omegacoi=(Pcoi/sum(d)+omegab)/omegab;
    else
        SEP=nan;
        omegacoi=nan;
    end
    % Test for calculation error
%     Pe_tmp=zeros(ngen,1);
%     for i=1:ngen
%         for j=1:ngen
%             ddelta=deltatt(i)-deltatt(j);
%             Pe_tmp(i)=Pe_tmp(i)+E(i)*E(j)*B(i,j)*sin(ddelta)+E(i)*E(j)*G(i,j)*cos(ddelta);
%         end
%     end
%     Pm-Pe_tmp-d*(omegacoi-1)*omegab

    
end
        
