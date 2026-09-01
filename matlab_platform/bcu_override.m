function out = bcu_override(arg)
% =========================================================================
% bcu_override —— 运行时参数覆盖通道（供 run_bcu_sweep 批量扫参使用）
%
% 用 persistent 保存一份「临时覆盖结构」。persistent 不受工作区 clear 影响
% （Cal_MM_CCT 的 clear 只清变量、不清 persistent），因此批量扫参时，
% 每组参数可以在这里设定，被 bcu_config() 合并进 cfg，一路注入到计算脚本，
% 全程不改任何文件。
%
% 用法：
%   bcu_override(struct('damping_ratio',[0.2;0.2;0.2]))  % 设定覆盖
%   ov = bcu_override();                                  % 读取当前覆盖([]表示无)
%   bcu_override('clear');                                % 清除覆盖
%
% 注意：单次入口 run_bcu 会在开头 bcu_override('clear')，确保单次运行只受
% bcu_config.m 文件控制；扫参结束后 run_bcu_sweep 也会清除，避免残留。
% =========================================================================
persistent OVER
if nargin == 0
    out = OVER;                      % 读取
    return;
end
if ischar(arg) || isstring(arg)
    if strcmp(char(arg),'clear'); OVER = []; end
    out = [];
elseif isstruct(arg)
    OVER = arg;                      % 设定
    out = [];
else
    error('bcu_override:badArg','参数须为 struct、''clear'' 或空。');
end
end
