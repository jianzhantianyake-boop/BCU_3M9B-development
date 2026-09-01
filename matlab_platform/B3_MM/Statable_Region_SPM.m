% =========================================================================
% 教学操作说明：结构保持模型的扩展平衡点与稳定区域搜索脚本
% 使用方法：
%  本文件是脚本。先在项目根目录运行 setup_bcu_paths()，再优先通过 run_bcu_beginner.m 选择对应模式。
% 参数：
%  没有函数形参；主要读取 base workspace 中的 preset、prefault、fault、postfault 等结构体，或读取本文件顶部的固定设置。
% 返回 / 工作区结果：
%  不返回函数值；计算结果保留在 base workspace，图形按原实现显示在 MATLAB 图窗中。
% 步骤：
%  1. 检查前置初始化与输入数据。 2. 执行本脚本的原始计算/搜索。 3. 保留结果变量或生成原始图窗。
% 单位：
%   角度通常为 rad，角速度为 rad/s，时间为 s，功率、电压和导纳通常为 pu；
%   个别中间变量为无量纲标志、迭代次数或矩阵索引，具体以变量定义为准。
% 前置条件：
%  需要保持原始脚本要求的 base workspace 状态；不要把它改写为函数，也不要跳过其上游初始化。
% 研究与验证边界：
%   本次只更新教学注释，不改变原始方程、参数、判据、求解器或绘图逻辑；
%   MATLAB 原生运行、收敛性和物理结论仍须在目标 MATLAB 环境中实际核验。
% =========================================================================
clear x_set x_set2 x_set3 ep_set 
system = "postfault";

N=2;
x=(-0.4:0.2:0.6)*2*pi;
n = length(x);
x_set = zeros(N,n^N);
for m = 0:(n^N - 1)
    for l = 1:N
        index = floor(m/n^(l-1));
        index = mod(index,n)+1;
        x_set(l,m+1) = x(index);
    end
end

x = (0:0.5:0.5)*2*pi;
n=length(x);
for m = 0:(n^6 - 1)
    for l = 1:6
        index = floor(m/n^(l-1));
        index = mod(index,n)+1;
        x_set2(l,m+1) = x(index);
    end
end

% y = (0.5:0.5:1);
% n=length(x);
% for m = 0:(n^6 - 1)
%     for l = 1:6
%         index = floor(m/n^(l-1));
%         index = mod(index,n)+1;
%         x_set3(l,m+1) = y(index);
%     end
% end
x_set3(:,1) = 0.5*ones(6,1);
x_set3(:,2) = 1*ones(6,1);


torralence = 1e-4; 
%% calculate EPs
m = 1;
ep_set = [];
options = optimoptions('fsolve','FunctionTolerance',1e-5,'display','off','MaxIterations',1000,'OptimalityTolerance',1e-5);
for n = 1:length(x_set(1,:))
    for p = 1:length(x_set2(1,:))
    for q = 1:length(x_set3(1,:))
    xep_extend = [x_set(:,n);x_set2(:,p);x_set3(:,q)];
    [xep_extend,ferr,exit,~,A] = fsolve(@f,xep_extend,options);

    xep = xep_extend(1:2);
    mmm=preset.m;
    xep23_1 = [xep(1)+(mmm(2)*xep(1)+mmm(3)*xep(2))/mmm(1); xep(2)+(mmm(2)*xep(1)+mmm(3)*xep(2))/mmm(1)];
    xep23_1 = mod(xep23_1,2*pi);
    delta_net_ep_1 = xep_extend(3:8)+(mmm(2)*xep(1)+mmm(3)*xep(2))/mmm(1);
    delta_net_ep_1 = mod(delta_net_ep_1,2*pi);
    coi_1 = (xep23_1(1)*mmm(2)+xep23_1(2)*mmm(3))/sum(mmm);
    xep = xep23_1-coi_1;
    delta_net_ep = delta_net_ep_1-coi_1;
    voltage_net_ep = xep_extend(9:14);

    if maxabs(ferr) < torralence && exit>0 && min(voltage_net_ep)>torralence
        if isnewxep(ep_set,xep,torralence)
           
            [V,Lambda]=eig(A);
            Lambda = diag(Lambda);
            sig = sign(sign(real(Lambda))+0.1); % zero counted as positive
            sig = (sig + 1)/2;                  % [0,1]
            flag = sum(sig);                    % number of non-negative eigenvalues

            v = V(:,~sig);                      % the stable sub-space
            
            ep_set(m).xep = xep; %#ok<*SAGROW> 
            ep_set(m).delta_net_ep = delta_net_ep; 
            ep_set(m).voltage_net_ep = voltage_net_ep; 
            ep_set(m).A = A;
            ep_set(m).Lambda = Lambda;
            ep_set(m).V = V;   
            ep_set(m).v = v;     % stable eigenvectors of unstable ep 
            ep_set(m).flag = flag;
            ep_set(m).xep_extend = [xep;delta_net_ep;voltage_net_ep];
           
            m = m+1;
            if flag == 0
            jacob=A;
            xeps=xep;
            end
        end
    end
    end
    end
end
%% Eigenvalue and Eigen vector
for n = 1:length(ep_set)
    A11 = ep_set(n).A(1:2,1:2);
    A12 = ep_set(n).A(1:2,3:14);
    A21 = ep_set(n).A(3:14,1:2);
    A22 = ep_set(n).A(3:14,3:14);
    A_reduce = A11-A12*inv(A22)*A21;
    [V,Lambda]=eig(A_reduce);
    Lambda = diag(Lambda);
    sig = sign(sign(real(Lambda))+0.1); % zero counted as positive
    sig = (sig + 1)/2;                  % [0,1]
    flag = sum(sig);                    % number of non-negative eigenvalues
    v = V(:,~sig);                      % the stable sub-space
    ep_set(n).flag_reduce = flag;
    ep_set(n).v_reduce = v;
    ep_set(n).V_reduce = V;
    ep_set(n).Lambda_reduce = Lambda;

    [V,Lambda]=eig(A22);
    Lambda = diag(Lambda);
    ep_set(n).Lambda_net = Lambda;
   

end


position1 = [1-mmm(2)/sum(mmm); -mmm(2)/sum(mmm)];
position2 = [-mmm(3)/sum(mmm); 1-mmm(3)/sum(mmm)];
position = [position1 position2];
position_net = [-mmm(2)/sum(mmm)*ones(6,1) -mmm(3)/sum(mmm)*ones(6,1)];
clear ep_set_ext
for n = 1:length(ep_set)
        m = (n-1)*6;
        ep_set_ext(m+1)=ep_set(n); %#ok<*AGROW> 
        ep_set_ext(m+2)=ep_set(n);
        ep_set_ext(m+3)=ep_set(n);
        ep_set_ext(m+4)=ep_set(n);
        ep_set_ext(m+5)=ep_set(n); %#ok<*AGROW> 
        ep_set_ext(m+6)=ep_set(n);
    
        ep_set_ext(m+2).xep = ep_set(n).xep - position*[2*pi;0];
        ep_set_ext(m+2).delta_net_ep = ep_set_ext(m+2).delta_net_ep - position_net*[2*pi;0];

        ep_set_ext(m+3).xep = ep_set(n).xep - position*[0;2*pi];
        ep_set_ext(m+3).delta_net_ep = ep_set_ext(m+3).delta_net_ep - position_net*[0;2*pi];

        ep_set_ext(m+4).xep = ep_set(n).xep - position*[2*pi;2*pi];
        ep_set_ext(m+4).delta_net_ep = ep_set_ext(m+4).delta_net_ep - position_net*[2*pi;2*pi];

        ep_set_ext(m+5).xep = ep_set(n).xep - position*[0;4*pi];
        ep_set_ext(m+5).delta_net_ep = ep_set_ext(m+5).delta_net_ep - position_net*[0;4*pi];

        ep_set_ext(m+6).xep = ep_set(n).xep - position*[2*pi;4*pi];
        ep_set_ext(m+6).delta_net_ep = ep_set_ext(m+6).delta_net_ep - position_net*[2*pi;4*pi];

end

for m = 1:length(ep_set)
    disp_v('Index',m);
    disp_v('Equilibrium',ep_set(m).xep);
    disp_v('Eigenvalue', ep_set(m).Lambda_reduce);
    disp_v('Eigenvector',ep_set(m).V_reduce);
end

%%
M = diag([ones(2,1); ones(12,1)*1e-15]);
options = odeset('Mass',M,'RelTol',1e-10,'AbsTol',[1e-8*ones(1,2),1e-12*ones(1,12)]);

figure;
hold on;
grid on;
color_code = {'blue','magenta','red','black'};
axis([-2*pi,2*pi,-2*pi,2*pi]);
for m = 1 : length(ep_set_ext)
        xep = ep_set_ext(m).xep;
        flag= ep_set_ext(m).flag_reduce;
        scatter(xep(1),xep(2),color_code{flag+1});
       
        if flag == 1
            v = ep_set_ext(m).v_reduce;
            perturb = 1e-2;
            perturb_xep_p = xep+v*perturb;
            perturb_xep_n = xep-v*perturb;
            theta_net = ep_set_ext(m).delta_net_ep;
            voltage_net = ep_set_ext(m).voltage_net_ep;
            temp_ini = [theta_net;voltage_net];
            delta1c= -mmm(2:ngen)'*perturb_xep_p/mmm(1);
            deltacc_p = [delta1c perturb_xep_p(1) perturb_xep_p(2)];
            delta1c= -mmm(2:ngen)'*perturb_xep_n/mmm(1);
            deltacc_n = [delta1c perturb_xep_n(1) perturb_xep_n(2)];
            [results_p,fval,exitflag,output,jacobian]=fsolve(@(x)Fun_AEfslove_SPM(x,deltacc_p,preset,"postfault"),temp_ini,optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off','TolX',1e-5));
            perturb_p = [perturb_xep_p;results_p];
            [results_n,fval,exitflag,output,jacobian]=fsolve(@(x)Fun_AEfslove_SPM(x,deltacc_n,preset,"postfault"),temp_ini,optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off','TolX',1e-5));
            perturb_n = [perturb_xep_n;results_n];

            [~ , x_p] = ode15s(@f_backward,[0,1],perturb_p,options);
            [~ , x_n] = ode15s(@f_backward,[0,1],perturb_n,options);
            x_all = [flip(x_n,1);x_p];
            plot(x_all(:,1),x_all(:,2),'k-','linewidth',1.5);
        end        
end


% %% trajectory
% X0=ep_set(2).xep+[0.4; 0];
% m=preset.m;
% delta1c= -m(2:ngen)'*X0(1:2)/m(1);
% deltacc = [delta1c X0(1) X0(2)];
% net_ini=[ep_set(2).delta_net_ep;ep_set(2).voltage_net_ep];
% [results,fval,exitflag,output,jacobian]=fsolve(@(x)Fun_AEfslove_SPM(x,deltacc,preset,"postfault"),net_ini,optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off','TolX',1e-5));
% xinit = [X0;net_ini-0.01];
% [tt , x_all] = ode15s(@f_backward,[0,3],xinit,options);
% plot(x_all(:,1),x_all(:,2),'k-','linewidth',1);





%%
function yes = isnewxep(ep_set,xep,torr)
    if isempty(ep_set)
        yes = 1;
        return;
    end
    minerr = inf;
    for m = 1 : length(ep_set)
        err = abs(xep - ep_set(m).xep);
        err = min(err, abs(2*pi-err));
        err = max(err);
        if minerr > err
            minerr = err;
        end
    end
    if(minerr>torr)
        yes = 1;
    else
        yes = 0;
    end
end
function disp_v(msg,v)
    disp([msg '=']);
    disp(v);
end
function dfdt = f(x)  
        dfdt = f_reducedstate_SPM(x); %dfdt = f_reducedstate2_SPM(x);
end
function dfdt = f_forward(t,x)
    dfdt = f_reducedstate_SPM(x);
end

function dfdt = f_backward(t,x)
    dfdt = f_reducedstate_SPM_backward(x);
end

function out = maxabs(in)

    out = abs(in);
    
    while length(out) > 1
        out = max(out);
    end

end

