function cfg = bcu_config()
% =========================================================================
% bcu_config —— BCU_3M9B 仿真平台的【唯一参数配置入口】
%
% 【怎么用】
%   1) 只修改本文件里 cfg.xxx = ... 右边的值，保存。
%   2) 回到命令窗口运行： run_bcu
%   3) run_bcu 会先校验配置、打印摘要、存参数快照，再跑你选的实验链路。
%
% 【为什么是函数而不是脚本变量】
%   原始 Cal_MM_CCT 会 clear 工作区。把配置写成函数，clear 清不掉它；
%   Cal_MM_Static / Cal_MM_CCT 在内部各自调用 bcu_config() 读取同一份真相，
%   因此参数在整条链路里保持一致，且不改变任何原始方程（见 override 钩子）。
%
% 【回归保证】
%   下面的默认值与原始 Cal_MM_Static.m / Cal_MM_CCT.m 完全一致；配置不改时，
%   结果与改造前 bit 级相同（用 verify/v0_baseline.m 可回归验证）。
%
% 单位约定：角度 rad，角速度 rad/s，时间 s，功率/电压/导纳 pu，频率 Hz。
% =========================================================================

% ======================== ① 实验选择 ====================================
% 决定调用哪条链路。可选值：
%   "reduced_cct"        网络约简模型：初始化 -> CCT/CUEP
%   "reduced_numerical"  网络约简模型：初始化 -> CCT/CUEP -> 数值轨迹
%   "reduced_region"     网络约简模型：初始化 -> 二维稳定区域
%   "spm_cct"            结构保持模型：初始化 -> CCT/CUEP
%   "spm_numerical"      结构保持模型：CCT/CUEP -> 数值轨迹
%   "spm_region"         结构保持模型：稳定区域
%   "two_machine_region_3d" / "two_machine_region_3d_gfl"  独立两机三维示例
% 注：参数注入当前完整覆盖 reduced_* 链路；SPM/两机链路暂用其脚本自带参数。
cfg.mode = "spm_cct";  % 选择实验链路

% ======================== ② 系统 / 案例 =================================
cfg.CaseName = 'case9_v2';   % MATPOWER case 文件名（不带 .m）。默认 3 机 9 节点。
cfg.f_base   = 60;           % 系统基频 Hz（omegab = 2*pi*f_base）。9-bus 用 60，39-bus 用 50。

% ======================== ③ 发电机参数（长度必须 == 发电机数）===========
% 顺序必须与 case 中发电机编号一致。默认对应 case9_v2 的 3 台机。
cfg.m             = [0.1254; 0.034; 0.016];   % 惯性系数 M=2H/omega_s (pu·s^2/rad)
cfg.damping_ratio = [0.1;    0.1;   0.1];     % 阻尼比 d_i/m_i；实际 d = m .* damping_ratio
cfg.Pm            = [0.8980; 1.3432; 0.9419]; % 机械输入功率 (pu)
cfg.xd1           = [0.0608; 0.1198; 0.1813]; % 暂态电抗 xd' (pu)
cfg.E             = [1.1083; 1.1071; 1.0606]; % 内电势幅值 (pu)
%   —— 切到 39-bus 时，把上面 5 个向量整组替换为 10 台机的值，并设
%      cfg.CaseName='case39_modified'; cfg.f_base=50; cfg.faultline 相应调整。

% ======================== ④ 故障设置 ====================================
cfg.faultline     = [9; 6];  % 故障支路 [FromBus; ToBus]，必须是 case 中存在的支路
cfg.faultposition = 0;       % 故障母线选择：0 取 faultline(1)，1 取 faultline(2)

% ======================== ⑤ 数值 / 求解器 ===============================
cfg.EquCal        = 2;       % 平衡点求解：1=自写 Newton，2=fsolve
cfg.PathEnergyCal = 0;       % 路径势能：0=Ray 近似，正整数 N=N 段梯形积分，-1=忽略
cfg.Tfault        = 2;       % 故障轨迹积分时长 (s)。需 > 预计 CCT。
cfg.Tunit         = 1e-4;    % 积分步长 (s)。越小越精确越慢（步长敏感性见 verify U5）。
cfg.TolFun        = 1e-50;   % fsolve 函数容差
cfg.TolX          = 1e-9;    % fsolve 步长容差

% ======================== ⑥ 验证 / 输出 =================================
cfg.run_selfcheck = true;    % 跑完 reduced_cct 后自动做 SEP/CUEP 残差 + LEA<=REA 自检
cfg.save_snapshot = true;    % 把本次配置存快照到 results/（可追溯，利于论文复现）

% ---- 运行时覆盖（批量扫参用）----
% run_bcu_sweep 通过 bcu_override(struct(...)) 设定的字段会在此覆盖上面的默认；
% 单次运行时 bcu_override() 为空，cfg 完全由本文件决定。
ov = bcu_override();
if isstruct(ov)
    fns = fieldnames(ov);
    for kk = 1:numel(fns)
        cfg.(fns{kk}) = ov.(fns{kk});
    end
end

end
