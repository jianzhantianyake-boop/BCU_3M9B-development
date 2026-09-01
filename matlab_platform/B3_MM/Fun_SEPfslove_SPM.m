% =========================================================================
% 教学操作说明：供 fsolve 求稳定平衡点使用的残差函数
% 使用方法：
%  按以下函数签名调用：function f= Fun_SEPfslove_SPM(x,preset,state,basevalue)
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
% x(2*nbus-ngen) = delta(1: ngen-1) | omegacoi(ngen)| delta_net(ngen+1:nbus) |
% V_net(nbus+1:2*nbus-ngen)
function f= Fun_SEPfslove_SPM(x,preset,state,basevalue)
% 教学分段提示：结构保持残差依次包含发电机角度/速度、网络角度、
% 网络电压和发电机/网络 P-Q 平衡；状态索引必须与 DAE 函数一致。
    Y_full=state.Yfull_mod;
    G=real(Y_full);    
    B=imag(Y_full);
    ngen=size(preset.genno,1);
    nbus=size(Y_full,1);
    m=preset.m;
    d=preset.d;
    Pm=preset.Pmpu;
    E=preset.Epu;
    omegab=basevalue.omegab;
    

    delta=zeros(ngen,1); % delta(ngen) is set as 0 as reference
    delta_net = zeros(nbus-ngen,1);
    V_net = zeros(nbus-ngen,1);
    Pe=zeros(ngen,1);
    Pnet = zeros(nbus-ngen,1);
    Qnet = zeros(nbus-ngen,1);
    domegacoi=x(ngen);

    % 第一组残差是 COI 角度、共同速度和网络状态的变量解包。
    for i=1:ngen-1
        delta(i)=x(i);
    end
    k=1;
    for i=ngen+1:nbus
        delta_net(k)=x(i);
        k=k+1;
    end
    k=1;
    for i=nbus+1:(2*nbus-ngen)
        V_net(k)=x(i);
        k=k+1;
    end
    clear k

    % Pe calculation
    % 第二组残差计算发电机有功平衡，第三/四组计算网络 P/Q 平衡。
    for i=1:ngen
        for j=1:ngen
            ddelta=delta(i)-delta(j);
            Pe(i)=Pe(i)+E(i)*E(j)*B(i,j)*sin(ddelta)+E(i)*E(j)*G(i,j)*cos(ddelta);
        end
        for l=1:(nbus-ngen)
            ddelta=delta(i)-delta_net(l);
            Pe(i)=Pe(i)+E(i)*V_net(l)*B(i,l+ngen)*sin(ddelta)+E(i)*V_net(l)*G(i,l+ngen)*cos(ddelta);
        end
    end

    PCOI=sum(Pm-Pe);   
    
    for i=1:ngen-1
        f(i)=Pm(i)-Pe(i)-m(i)/sum(m)*PCOI+m(i)/sum(m)*sum(d)*domegacoi-d(i)*domegacoi;
    end
    f(ngen)= sum(Pm-Pe)-sum(d)*domegacoi;    % the nth equ

    % P calculation of Bus
    for i=1:(nbus-ngen)
        for j=1:ngen
            ddelta=delta_net(i)-delta(j);
            Pnet(i)=Pnet(i)+V_net(i)*E(j)*B(i+ngen,j)*sin(ddelta)+V_net(i)*E(j)*G(i+ngen,j)*cos(ddelta);
        end
        for l=1:(nbus-ngen)
            ddelta=delta_net(i)-delta_net(l);
            Pnet(i)=Pnet(i)+V_net(i)*V_net(l)*B(i+ngen,l+ngen)*sin(ddelta)+V_net(i)*V_net(l)*G(i+ngen,l+ngen)*cos(ddelta);
        end
        for h=1:size(preset.Sload,1)
            if (preset.Sload(h,1)==state.Transform(i+ngen))
                Pnet(i)=Pnet(i)+preset.Sload(h,2);
            end
        end
    end

    for i=1:(nbus-ngen)
        f(ngen+i)=Pnet(i);
    end

    % Q calculation of Bus
    for i=1:(nbus-ngen)
        for j=1:ngen
            ddelta=delta_net(i)-delta(j);
            Qnet(i)=Qnet(i)-E(j)*B(i+ngen,j)*cos(ddelta)+E(j)*G(i+ngen,j)*sin(ddelta);
        end
        for l=1:(nbus-ngen)
            ddelta=delta_net(i)-delta_net(l);
            Qnet(i)=Qnet(i)-V_net(l)*B(i+ngen,l+ngen)*cos(ddelta)+V_net(l)*G(i+ngen,l+ngen)*sin(ddelta);
        end
        for h=1:size(preset.Sload,1)
            if (preset.Sload(h,1)==state.Transform(i+ngen))
                Qnet(i)=Qnet(i)+preset.Sload(h,3)/V_net(i);
            end
        end
    end

    for i=1:(nbus-ngen)
        f(nbus+i)=Qnet(i);
    end

end
