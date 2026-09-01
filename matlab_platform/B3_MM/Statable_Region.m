% =========================================================================
% 教学操作说明：网络约简模型的二维平衡点与稳定区域搜索脚本
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
N=2;
x=(-1:0.1:1)*2*pi;
n = length(x);
x_set = zeros(N,n^N);
for m = 0:(n^N - 1)
    for l = 1:N
        index = floor(m/n^(l-1));
        index = mod(index,n)+1;
        x_set(l,m+1) = x(index);
    end
end
torralence = 1e-2; 
clear ep_set ep_set_ext
%% calculate EPs
m = 1;
ep_set = [];
options = optimoptions('fsolve','FunctionTolerance',1e-10,'MaxIterations',10000,'OptimalityTolerance',1e-10);
for n = 1:length(x_set(1,:))
    xep = x_set(:,n);
    [xep,ferr,~,~,A] = fsolve(@f,xep,options);

    mmm=preset.m;
    xep23_1 = [xep(1)+(mmm(2)*xep(1)+mmm(3)*xep(2))/mmm(1); xep(2)+(mmm(2)*xep(1)+mmm(3)*xep(2))/mmm(1)];
    xep23_1 = mod(xep23_1,2*pi);
    coi_1 = (xep23_1(1)*mmm(2)+xep23_1(2)*mmm(3))/sum(mmm);
    xep = xep23_1-coi_1;
    
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

position1 = [1-mmm(2)/sum(mmm); -mmm(2)/sum(mmm)];
position2 = [-mmm(3)/sum(mmm); 1-mmm(3)/sum(mmm)];
position = [position1 position2];

for n = 1:length(ep_set)
        m = (n-1)*6;
        ep_set_ext(m+1)=ep_set(n); %#ok<*AGROW> 
        ep_set_ext(m+2)=ep_set(n);
        ep_set_ext(m+3)=ep_set(n);
        ep_set_ext(m+4)=ep_set(n);
        ep_set_ext(m+5)=ep_set(n); %#ok<*AGROW> 
        ep_set_ext(m+6)=ep_set(n);
    
        ep_set_ext(m+2).xep = ep_set(n).xep - position*[2*pi;0];

    
        ep_set_ext(m+3).xep = ep_set(n).xep - position*[0;2*pi];
    
        ep_set_ext(m+4).xep = ep_set(n).xep - position*[2*pi;2*pi];

        ep_set_ext(m+5).xep = ep_set(n).xep - position*[0;4*pi];

        ep_set_ext(m+6).xep = ep_set(n).xep - position*[2*pi;4*pi];

end
f1=figure(1);
figure(f1);
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
%% 
% xep = ep_set_ext(19).xep;
% v = ep_set_ext(19).v;
% perturb = 1e-2;
% [~ , x_p] = ode45(@f_backward,[0,50],xep+v*perturb,odeset('RelTol',1e-5));
% [~ , x_n] = ode45(@f_backward,[0,50],xep-v*perturb,odeset('RelTol',1e-5));
% x_all = [flip(x_n,1);x_p];
% plot(x_all(:,1),x_all(:,2),'b-','linewidth',1);
% axis([0,2.5,0,3.5]);
%%

for m = 1:length(ep_set)
    disp_v('Index',m);
    disp_v('Equilibrium',ep_set(m).xep);
    disp_v('Eigenvalue', ep_set(m).Lambda);
    disp_v('Eigenvector',ep_set(m).V);
end

function yes = isnewxep(ep_set,xep,torr)
    if isempty(ep_set)
        yes = 1;
        return;
    end
    minerr = inf;
    for m = 1 : length(ep_set)
        err = abs(xep - ep_set(m).xep);
        %err = min(err, abs(2*pi-err));
        %err = max(err);
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
        dfdt = f_reducedstate(x);
end
function dfdt = f_forward(t,x)
    dfdt = f(x);
end

function dfdt = f_backward(t,x)
    dfdt = -f(x);
end
function out = maxabs(in)

    out = abs(in);
    
    while length(out) > 1
        out = max(out);
    end

end
