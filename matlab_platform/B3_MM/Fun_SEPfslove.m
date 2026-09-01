% =========================================================================
% 教学操作说明：供 fsolve 求稳定平衡点使用的残差函数
% 使用方法：
%  按以下函数签名调用：function f= Fun_SEPfslove(delta_omega,preset,state,basevalue)
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
%% unit: omega in rad/s (in SRF)
function f= Fun_SEPfslove(delta_omega,preset,state,basevalue)
% 教学分段提示：delta_omega 的最后一项是共同速度，前面的角度采用
% 相对最后一台发电机的表示；末尾再用 COI 约束恢复全机角度。
    Y_red=state.Yred;
    G=real(Y_red);    
    B=imag(Y_red);
    ngen=size(Y_red,1);
    m=preset.m;
    d=preset.d;
    Pm=preset.Pmpu;
    E=preset.Epu;
    omegab=basevalue.omegab;
    

% f=State.Pm*0;
    delta=zeros(ngen,1);
    Pe=zeros(ngen,1);
        omega=delta_omega(ngen);
        for i=1:ngen-1
            delta(i)=delta_omega(i);
        end

    % Pe calculation
        for i=1:ngen
            for j=1:ngen
                ddelta=delta(i,1)-delta(j,1);
                Pe(i)=Pe(i)+E(i)*E(j)*B(i,j)*sin(ddelta)+E(i)*E(j)*G(i,j)*cos(ddelta);
            end
        end

        PCOI=sum(Pm-Pe)-sum(d)*omega;

    % 先构造 fsolve 的角度和速度残差，再计算 Yred 上的电功率失配。
    for i=1:ngen-1
        f(i)=Pm(i)-Pe(i)-m(i)/sum(m)*PCOI-d(i)*omega;
    end
    f(ngen)= sum(Pm-Pe)-sum(d)*omega;    % the nth equ

end
