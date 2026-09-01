% =========================================================================
% 教学操作说明：两机示例模型的状态导数计算函数
% 使用方法：
%  按以下函数签名调用：function dfdt = f_2m(x)
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
function dfdt = f_2m(x)


delta12 = x(1);
omega12 = x(2);
omegasum = x(3);

G1 = evalin('base','G1');
G2 = evalin('base','G2');
G12 = evalin('base','G12');
B12 = evalin('base','B12');
E1 = evalin('base','E1');
E2 = evalin('base','E2');
Pm1 = evalin('base','Pm1');
Pm2 = evalin('base','Pm2');
H1 = evalin('base','H1');
H2 = evalin('base','H2');
D1 = evalin('base','D1');
D2 = evalin('base','D2');

Pe1 = E1^2*(G12+G1) + E1*E2*G12*cos(delta12) + E1*E2*B12*sin(delta12);
Pe2 = E2^2*(G12+G2) + E1*E2*G12*cos(delta12) - E1*E2*B12*sin(delta12);

dfdt(1) = omega12;% delta12
dfdt(2) = Pm1/H1 - Pm2/H2 - Pe1/H1 + Pe2/H2 - omega12*(D1/H1+D2/H2)/2 - omegasum*(D1/H1-D2/H2)/2;%omega12
dfdt(3) = Pm1/H1 + Pm2/H2 - Pe1/H1 - Pe2/H2 - omega12*(D1/H1-D2/H2)/2 - omegasum*(D1/H1+D2/H2)/2;%omegasum


dfdt = dfdt.';

end
