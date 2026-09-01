% =========================================================================
% 教学操作说明：从既有 MAT 数据文件读取初值并绘制三维状态示意图的辅助脚本
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
close all
clear
%% Data Import
load('Data_InitState_Dp20_H0.1_ZL1.mat');

%% Critical Stable Point Identification
Stb.up.state=[];
Stb.down.state=[];
D1.pos.state=[];
D1.neg.state=[];
D2.pos.state=[];
D2.neg.state=[];
for k3=1:n3
        thetapg_tmp=thetapg0(1,1,k3);
    for k1=1:n1
        omegapg_tmp=omegapg0(k1,1,1);
        Stb.down.flag=0;
        Stb.up.flag=0;
        for k2=1:n2
        %% Stable Boundary
            omegagdown_tmp=omegag0(1,k2,1);
            omegagup_tmp=omegag0(1,n2+1-k2,1);
            % search for downward omegag0
            if(flag_unstb(k1,k2,k3)==0&&Stb.down.flag==0)
                statedown_tmp=[thetapg_tmp,omegapg_tmp,omegagdown_tmp];
                Stb.down.flag=1;
                n_stbdown=k2;
                if(isempty(Stb.down.state)==1)
                    Stb.down.state=statedown_tmp;
                else
                    Stb.down.state=[Stb.down.state;statedown_tmp];
                end
            end
            % search for downward omegag0
            if(flag_unstb(k1,n2+1-k2,k3)==0&&Stb.up.flag==0)
                stateup_tmp=[thetapg_tmp,omegapg_tmp,omegagup_tmp];
                Stb.up.flag=1;
                n_stbup=n2+1-k2;
                if(isempty(Stb.up.state)==1)
                    Stb.up.state=stateup_tmp;
                else
                    Stb.up.state=[Stb.up.state;stateup_tmp];
                end
            end
        %% Wd1 Boudary
        omegag_tmp=omegag0(1,k2,1);
        if(flag_d1(k1,k2,k3)==1&&flag_unstb(k1,k2,k3)==0)    % stable boundary is found
                stated1p_tmp=[thetapg_tmp,omegapg_tmp,omegag_tmp];
                if(isempty(D1.pos.state)==1)
                    D1.pos.state=stated1p_tmp;
                else
                    D1.pos.state=[D1.pos.state;stated1p_tmp];
                end
        elseif(flag_d1(k1,k2,k3)==2&&flag_unstb(k1,k2,k3)==0)
                stated1n_tmp=[thetapg_tmp,omegapg_tmp,omegag_tmp];
                if(isempty(D1.neg.state)==1)
                    D1.neg.state=stated1n_tmp;
                else
                    D1.neg.state=[D1.neg.state;stated1n_tmp];
                end
        end
        %% Wd2 Boudary
        if(flag_d2(k1,k2,k3)==1&&flag_unstb(k1,k2,k3)==0)    % stable boundary is found
                stated2p_tmp=[thetapg_tmp,omegapg_tmp,omegag_tmp];
                if(isempty(D2.pos.state)==1)
                    D2.pos.state=stated2p_tmp;
                else
                    D2.pos.state=[D2.pos.state;stated2p_tmp];
                end
        elseif(flag_d2(k1,k2,k3)==2&&flag_unstb(k1,k2,k3)==0)
                stated2n_tmp=[thetapg_tmp,omegapg_tmp,omegag_tmp];
                if(isempty(D2.neg.state)==1)
                    D2.neg.state=stated2n_tmp;
                else
                    D2.neg.state=[D2.neg.state;stated2n_tmp];
                end
        end
        end
    end
end

%% Plot
figure(1);
    x_down=Stb.down.state(:,1);
    y_down=Stb.down.state(:,2);
    z_down=Stb.down.state(:,3);
    x_up=Stb.up.state(:,1);
    y_up=Stb.up.state(:,2);
    z_up=Stb.up.state(:,3);
% point
%     scatter3(x_down,y_down,z_down,'MarkerEdgeColor','c'); hold on;
%     scatter3(x_up,y_up,z_up,'MarkerEdgeColor','m'); hold on;
% surface
    [X_down,Y_down,Z_down]=griddata(x_down,y_down,z_down,linspace(min(x_down),max(x_down))',linspace(min(y_down),max(y_down)),'v4');
    surf(X_down,Y_down,Z_down,'FaceAlpha',0.1,'EdgeColor','none','FaceColor','y'); hold on;
    [X_up,Y_up,Z_up]=griddata(x_up,y_up,z_up,linspace(min(x_up),max(x_up))',linspace(min(y_up),max(y_up)),'v4');
    surf(X_up,Y_up,Z_up,'FaceAlpha',0.1,'EdgeColor','none','FaceColor','y');
%% D1 （GFL damping）
%     x_d1p=D1.pos.state(:,1);
%     y_d1p=D1.pos.state(:,2);
%     z_d1p=D1.pos.state(:,3);
%     x_d1n=D1.neg.state(:,1);
%     y_d1n=D1.neg.state(:,2);
%     z_d1n=D1.neg.state(:,3);
%     scatter3(x_d1p,y_d1p,z_d1p,5,'MarkerEdgeColor','c','Marker','*');  hold on;
%     scatter3(x_d1n,y_d1n,z_d1n,5,'MarkerEdgeColor','m','Marker','*');  hold on;
%% D2 (GFM damping)
    x_d2p=D2.pos.state(:,1);
    y_d2p=D2.pos.state(:,2);
    z_d2p=D2.pos.state(:,3);
    x_d2n=D2.neg.state(:,1);
    y_d2n=D2.neg.state(:,2);
    z_d2n=D2.neg.state(:,3);
    scatter3(x_d2p,y_d2p,z_d2p,5,'MarkerEdgeColor','c','Marker','o');  hold on;
    scatter3(x_d2n,y_d2n,z_d2n,5,'MarkerEdgeColor','m','Marker','o');  hold on;
%% 2D--omegag=1 plane
    deltamin=thetapg0(1,1,1);
    deltamax=thetapg0(1,1,n3);
    omegapgmin=omegapg0(1,1,1);
    omegapgmax=omegapg0(n1,1,1);
    [x,y]=meshgrid(deltamin:(deltamax-deltamin)/40:deltamax,omegapgmin:(omegapgmax-omegapgmin)/40:omegapgmax);
    z=ones(41,41);
    surf(x,y,z,'FaceAlpha',0.1,'EdgeColor','none','FaceColor','b');
    clear x y z
