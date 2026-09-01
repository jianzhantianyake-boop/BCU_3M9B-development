% =========================================================================
% 教学操作说明：独立两机参数下的三维平衡点与稳定区域探索脚本
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
clear;

%% Parameter
Z1 = 0+0.5j;
Z2 = 0+0.3j;
Zl = 0.2844+0.0306j;
Pm1 = 2.2;
Pm2 = 1.3;
H1 = 0.6;
H2 = 0.5;
D1 = H1*1;
D2 = H2*0.9;
Y12 = 1/(Z1+Z2+Z1*Z2/Zl);
Y1 = 1/(Z1+Zl+Z1*Zl/Z2);
Y2 = 1/(Z2+Zl+Z2*Zl/Z1);
G1 = real(Y1);
G2 = real(Y2);
G12 = -real(Y12);
B12 = -imag(Y12);
G1_fault = real(1/Z1);
G2_fault = real(1/Z2);
E1 = 1; E2 = 1;

%%
global N;
N=3;
try
    N;
catch
    N = 2;
end
x=(0:0.1:1)*2*pi;
n = length(x);
x_set = zeros(N,n);
for mm = 1:n
   x_set(1,mm) = x(mm);
end
torralence = 1e-2; 
%% calculate EPs
m = 1;
ep_set = [];
options = optimoptions('fsolve','FunctionTolerance',1e-10,'MaxIterations',10000,'OptimalityTolerance',1e-10);
for n = 1:length(x_set(1,:))
    xep = x_set(:,n);
    [xep,ferr,~,~,A] = fsolve(@f1,xep,options);

    %xep = mod(xep,2*pi);
    
    if maxabs(ferr) < torralence
        if isnewxep(ep_set,xep,torralence)
           
            [V,Lambda]=eig(A);
            Lambda = diag(Lambda);
            sig = sign(sign(real(Lambda))+0.1); % zero counted as positive
            sig = (sig + 1)/2;                  % [0,1]
            flag = sum(sig);                    % number of non-negative eigenvalues

            v = V(:,~sig);                      % the stable sub-space
            
            ep_set(m).xep = real(xep); %#ok<*SAGROW> 
            ep_set(m).A = A;
            ep_set(m).Lambda = Lambda;
            ep_set(m).V = V;   
            ep_set(m).v = v;     % stable eigenvectors of unstable ep 
            ep_set(m).flag = flag;
           
            m = m+1;
            if flag == 0
            jacob=A;
            xeps=xep;
            end
        end
    end
end

for m = 1:length(ep_set)
    disp_v('Index',m);
    disp_v('Equilibrium',ep_set(m).xep);
    disp_v('Eigenvalue', ep_set(m).Lambda);
    disp_v('Eigenvector',ep_set(m).V);
end
if N == 2
for n = 1:length(ep_set)
        m = (n-1)*3;
        ep_set_ext(m+1)=ep_set(n); %#ok<*AGROW> 
        ep_set_ext(m+2)=ep_set(n);
        ep_set_ext(m+3)=ep_set(n);
    
        ep_set_ext(m+2).xep(1) = ep_set(n).xep(1) - 2*pi;
    
        ep_set_ext(m+3).xep(1) = ep_set(n).xep(1) + 2*pi;
    
end
figure;
hold on;
grid on;
color_code = {'blue','magenta','red','black'};
axis([-2*pi,2*pi,-2*pi,2*pi]);
for m = 1 : length(ep_set_ext)
        xep = ep_set_ext(m).xep;
        flag= ep_set_ext(m).flag;
        scatter(xep(1),xep(2),color_code{flag+1});
       
        if flag == 1
            v = ep_set_ext(m).v;
            perturb = 1e-2;
            [~ , x_p] = ode45(@f_backward,[0,50],xep+v*perturb,odeset('RelTol',1e-5));
            [~ , x_n] = ode45(@f_backward,[0,50],xep-v*perturb,odeset('RelTol',1e-5));
            x_all = [flip(x_n,1);x_p];
            plot(x_all(:,1),x_all(:,2),'k-','linewidth',1.5);
        end        
end
elseif N == 3
for n = 1:length(ep_set)
        m = (n-1)*3;
        ep_set_ext(m+1)=ep_set(n); %#ok<*AGROW> 
        ep_set_ext(m+2)=ep_set(n);
        ep_set_ext(m+3)=ep_set(n);
        ep_set_ext(m+2).xep(1) = ep_set(n).xep(1) - 2*pi;
        ep_set_ext(m+3).xep(1) = ep_set(n).xep(1) + 2*pi;   
end
figure;

color_code = {'blue','magenta','green','red','cyan','blue','magenta','green','red','cyan'};
color_code2 = {[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],[0.8500 0.3250 0.0980],'red','cyan','blue','magenta','green','red','cyan'};
for m = 1:length(ep_set_ext)
    xep = ep_set_ext(m).xep;
    flag = ep_set_ext(m).flag;
    scatter3(xep(1),xep(2),xep(3),color_code{flag+1});
%     xlabel('delta12');
%     ylabel('omega12');
%     zlabel('omega_{sum}');
    hold on;
end    
n = 1;
for m = 1 : length(ep_set_ext)        
    flag = ep_set_ext(m).flag;
    if flag 
        xep = ep_set_ext(m).xep;
        v = ep_set_ext(m).v;
        perturb = 5e-1;

        if flag == 1
            for alpha = (0:0.005:1)*2*pi
                vp = v(:,1)*sin(alpha) + v(:,2)*cos(alpha);
                [~ , x_all] = ode45(@f_backward,[0,5],xep+vp*perturb);
                plot3(x_all(:,1),x_all(:,2),x_all(:,3),'color',color_code2{n});
            end
            n = n + 1;
        elseif flag == 2
            for beta = [-1,1]
                vp = v*beta;
                [~ , x_all] = ode45(@f_backward,[0,1],xep+vp*perturb);
                plot3(x_all(:,1),x_all(:,2),x_all(:,3),'color',color_code2{n});
            end
            n = n + 1;
        end
        
    end
end
% % plot trajectory
 [tt , x_all] = ode78(@f_forward,[0,200],[-10, 3, -20],odeset('RelTol',1e-5));
% [tt2 , x_all2] = ode78(@f_forward,[0,200],[0.311-4*pi, 4.75, 10],odeset('RelTol',1e-5));
% [tt3 , x_all3] = ode78(@f_forward,[0,200],[-1.16, -9.67, -20],odeset('RelTol',1e-5));
% [tt4 , x_all4] = ode78(@f_forward,[0,200],[1.02, 2.11, 10],odeset('RelTol',1e-5));
% [tt5 , x_all5] = ode78(@f_forward,[0,200],[1.02, 1.11, -20],odeset('RelTol',1e-5));
% [tt6 , x_all6] = ode78(@f_forward,[0,1200],[-2.667, 0.045+1, -60],odeset('RelTol',1e-5));%[-2.667, 0.045+1, 80]
% [tt7 , x_all7] = ode78(@f_forward,[0,200],[-1, -1.807, 20],odeset('RelTol',1e-5));
% [tt8 , x_all8] = ode78(@f_forward,[0,200],[3, 0, 0],odeset('RelTol',1e-5));

plot3(x_all(:,1),x_all(:,2),x_all(:,3),'black','linewidth',1.5);
% plot3(x_all2(:,1),x_all2(:,2),x_all2(:,3),'black','linewidth',1.5);
% plot3(x_all3(:,1),x_all3(:,2),x_all3(:,3),'black','linewidth',1.5);
% plot3(x_all4(:,1),x_all4(:,2),x_all4(:,3),'black','linewidth',1.5);
% plot3(x_all5(:,1),x_all5(:,2),x_all5(:,3),'black','linewidth',1.5);
% plot3(x_all6(:,1),x_all6(:,2),x_all6(:,3),'black','linewidth',1.5);
% plot3(x_all7(:,1),x_all7(:,2),x_all7(:,3),'black','linewidth',1.5);
% plot3(x_all8(:,1),x_all8(:,2),x_all8(:,3),'black','linewidth',1.5);

% plot trajectory
t_fault=1.65; %1.3;
[tt , x_all] = ode78(@f_fault,[0,t_fault],ep_set_ext(1).xep',odeset('RelTol',1e-5));
[tt2 , x_all2] = ode78(@f_forward,[t_fault,40],[x_all(end,1), x_all(end,2),x_all(end,3)],odeset('RelTol',1e-5));
plot3(x_all(:,1),x_all(:,2),x_all(:,3),'red','linewidth',2);
plot3(x_all2(:,1),x_all2(:,2),x_all2(:,3),'blue','linewidth',2);




end

%%
H1 = 0.6;
H2 = 0.5;
D1 = H1*1;
D2 = H2*1;

m = 1;
ep_set = [];
options = optimoptions('fsolve','FunctionTolerance',1e-10,'MaxIterations',10000,'OptimalityTolerance',1e-10);
for n = 1:length(x_set(1,:))
    xep = x_set(:,n);
    [xep,ferr,~,~,A] = fsolve(@f1,xep,options);

    %xep = mod(xep,2*pi);
    
    if maxabs(ferr) < torralence
        if isnewxep(ep_set,xep,torralence)
           
            [V,Lambda]=eig(A);
            Lambda = diag(Lambda);
            sig = sign(sign(real(Lambda))+0.1); % zero counted as positive
            sig = (sig + 1)/2;                  % [0,1]
            flag = sum(sig);                    % number of non-negative eigenvalues

            v = V(:,~sig);                      % the stable sub-space
            
            ep_set(m).xep = real(xep); %#ok<*SAGROW> 
            ep_set(m).A = A;
            ep_set(m).Lambda = Lambda;
            ep_set(m).V = V;   
            ep_set(m).v = v;     % stable eigenvectors of unstable ep 
            ep_set(m).flag = flag;
           
            m = m+1;
            if flag == 0
            jacob=A;
            xeps=xep;
            end
        end
    end
end


for n = 1:length(ep_set)
        m = (n-1)*3;
        ep_set_ext(m+1)=ep_set(n); %#ok<*AGROW> 
        ep_set_ext(m+2)=ep_set(n);
        ep_set_ext(m+3)=ep_set(n);
        ep_set_ext(m+2).xep(1) = ep_set(n).xep(1) - 2*pi;
        ep_set_ext(m+3).xep(1) = ep_set(n).xep(1) + 2*pi;   
end

color_code = {'blue','magenta','green','red','cyan','blue','magenta','green','red','cyan'};
color_code2 = {[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],[0.8500 0.3250 0.0980],'red','cyan','blue','magenta','green','red','cyan'};
for m = 1:length(ep_set_ext)
    xep = ep_set_ext(m).xep;
    flag = ep_set_ext(m).flag;
%     scatter3(xep(1),xep(2),xep(3),color_code{flag+1});
    hold on;
end    
n = 1;
for m = 1 : length(ep_set_ext)        
    flag = ep_set_ext(m).flag;
    if flag 
        xep = ep_set_ext(m).xep;
        v = ep_set_ext(m).v;
        perturb = 5e-1;

        if flag == 1
            for alpha = (0:0.005:1)*2*pi
                vp = v(:,1)*sin(alpha) + v(:,2)*cos(alpha);
                [~ , x_all] = ode45(@f_backward,[0,10],xep+vp*perturb);
                plot3(x_all(:,1),x_all(:,2),x_all(:,3),'k');
            end
            n = n + 1;
        elseif flag == 2
            for beta = [-1,1]
                vp = v*beta;
                [~ , x_all] = ode45(@f_backward,[0,2],xep+vp*perturb);
                plot3(x_all(:,1),x_all(:,2),x_all(:,3),'k');
            end
            n = n + 1;
        end
        
    end
end

axis([-2*pi 2*pi -15 10 -100 100]);
xticks(-2*pi:pi/2:2*pi);
xticklabels({'$-2\pi$', '', '$-\pi$', '','$0$', '','$\pi$', '','$2\pi$'});
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 14);
xlabel('\delta_{12}');
ylabel('\omega_{12}');
zlabel('\omega_{sum}');





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
function dfdt = f1(x)  
 global N; 
    switch N
        case 2
            dfdt = f_2m_reduce(x);
        case 3
            dfdt = f_2m(x);
    end 
end
function dfdt = f_forward(t,x)
    dfdt = f1(x);
end

function dfdt = f_backward(t,x)
    dfdt = -f1(x);
end

function dfdt = f_fault(t,x)
 global N; 
    switch N
        case 2
            dfdt = f_2m_reduce_fault(x);
        case 3
            dfdt = f_2m_fault(x);
    end 
end

function out = maxabs(in)
    out = abs(in);
    while length(out) > 1
        out = max(out);
    end
end
