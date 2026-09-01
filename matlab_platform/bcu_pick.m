function v = bcu_pick(cfg, field, default)
% =========================================================================
% bcu_pick —— 参数覆盖钩子（override hook）
%
% 若配置结构体 cfg 提供了非空字段 field，则返回该值（用户覆盖）；
% 否则返回 default（原始默认）。这是让 bcu_config 覆盖原始脚本参数、
% 又保证「配置缺字段时仍等价于原始默认」的唯一机制。
%
% 用法（在 Cal_MM_Static / Cal_MM_CCT 内部）：
%   preset.m = bcu_pick(BCUCFG, 'm', [0.1254;0.034;0.016]);
% 第三个参数 default 保留了原始默认值，因此本改造对原代码零回归。
% =========================================================================
if isstruct(cfg) && isfield(cfg, field) && ~isempty(cfg.(field))
    v = cfg.(field);
else
    v = default;
end
end
