% =========================================================================
% 教学操作说明：读取已存在 Group 结构体并叠加稳定边界与 CUEP 曲线的辅助绘图脚本
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
for i=1:10
        figure(1);
        plot(Group.delta_stb(1:5e4,i),'Color',[(50+20*i)/255 150/255 (250-20*i)/255]); hold on;
%         strlb(no_group)="group"+no_group;
%         legend(strlb);
        cycle=size(Group.delta_stb(1:5e4,i),1);
        plot(Group.CUEP(i)*ones(cycle,1),'Color',[(50+20*i)/255 150/255 (250-20*i)/255],'LineStyle',':');    hold on;
%         strUEP(no_group)="CUEP"+no_group;
%         legend(strUEP);
%         figure(10+no_group);
%         plot(Group.omega_stb(:,no_group)); hold on;
%         figure(3)
%         plot(Group.delta_unstb(:,no_group)); hold on;
%         strlb(no_group)="group"+no_group;
%         legend(strlb);
%         plot(Group.CUEP(no_group)*ones(cycle_stb,1),':');    hold on;
%         strUEP(no_group)="CUEP"+no_group;
%         legend(strUEP);
    end
