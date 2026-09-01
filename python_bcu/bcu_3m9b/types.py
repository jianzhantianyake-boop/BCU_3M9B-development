"""平台数据结构说明书。

使用方法:
    先用 ``CaseData`` 装入案例数据,调用潮流得到 ``PowerFlowResult``,再由
    ``to_pfdata`` 整理成 ``PFData``; 随后用 ``Preset``、``BaseValue`` 配置动态
    参数与基值,网络处理产出 ``NetworkState``,仿真产出 ``Trajectory``,一次
    静态初始化的全部结果打包进 ``StaticResult``,CCT 结果打包进 ``CCTResult``。

设计说明:
    MATLAB 版本大量依赖结构体和 ``evalin('base')``。这里改用 dataclass 把输入、
    输出和中间状态显式写出,便于在实验中保存、检查和复用,不再隐式读全局变量。
"""

from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import numpy as np


@dataclass
class CaseData:
    """MATPOWER 案例的最小数据容器。

    使用方法:
        直接构造并传入下面各字段; ``bus``、``gen``、``branch`` 均沿用 MATPOWER
        的原始列顺序,母线编号从 1 开始计数。
    参数:
        base_mva:系统功率基值(MVA)。
        bus:母线数据矩阵。
        gen:发电机数据矩阵。
        branch:支路数据矩阵。
        gencost:可选的发电成本矩阵。
        name:案例名称,用于打印识别。
    """

    base_mva: float
    bus: np.ndarray
    gen: np.ndarray
    branch: np.ndarray
    gencost: Optional[np.ndarray] = None
    name: str = "Unnamed case"


@dataclass
class PowerFlowResult:
    """纯 NumPy 交流牛顿潮流的返回容器。

    使用方法:
        读取 ``solve_power_flow`` 的返回值即可; ``success`` 判断是否收敛,
        ``residual_norm`` 查看最终功率残差范数。
    参数:
        case:本次求解使用的 ``CaseData``。
        voltage:各母线复电压(pu)。
        bus_injection_mva:各母线复功率注入(MVA)。
        gen:反算后的发电机出力矩阵。
        branch_flow:支路两端潮流。
        success:是否收敛。
        iterations:牛顿迭代次数。
        residual_norm:最终残差二范数。
    """

    case: CaseData
    voltage: np.ndarray
    bus_injection_mva: np.ndarray
    gen: np.ndarray
    branch_flow: np.ndarray
    success: bool
    iterations: int
    residual_norm: float


@dataclass
class PFData:
    """对应 MATLAB ``Fun_ResultBack`` 输出的项目内部潮流结构。

    使用方法:
        由 ``to_pfdata`` 生成,供网络约简和结构保持模型读取; 通过 ``nbus``、
        ``ngen``、``nload`` 属性直接取规模。
    参数说明:
        sbase 为功率基值; voltage 为 [幅值, 角度] 两列; bus_pq 为母线 P/Q; 
        gen_no/gen_pq/gen_voltage 为发电机编号、出力和端电压; load_* 同理为负荷; 
        rxb 为线路 [i, j, r, x, b]; branch_powerflow 为支路潮流。
    """

    sbase: float
    voltage: np.ndarray
    bus_pq: np.ndarray
    gen_no: np.ndarray
    gen_pq: np.ndarray
    gen_voltage: np.ndarray
    load_no: np.ndarray
    load_pq: np.ndarray
    load_voltage: np.ndarray
    rxb: np.ndarray
    branch_powerflow: np.ndarray

    @property
    def nbus(self) -> int:
        # 用法:取母线总数(voltage 的行数)。
        return int(self.voltage.shape[0])

    @property
    def ngen(self) -> int:
        # 用法:取发电机台数(gen_no 的元素个数)。
        return int(self.gen_no.size)

    @property
    def nload(self) -> int:
        # 用法:取负荷母线个数(load_no 的元素个数)。
        return int(self.load_no.size)


@dataclass
class BaseValue:
    """系统基值容器。

    使用方法:
        用 ``BaseValue(sbase=...)`` 覆盖默认功率基值; ``omega_b`` 为基准角速度
        (rad/s),默认按 60 Hz 计算。
    """

    omega_b: float = 2.0 * np.pi * 60.0
    sbase: float = 100.0


@dataclass
class Preset:
    """三机九母线动态参数容器。

    使用方法:
        构造时传入 m、d、pmpu、xd1、epu,其余字段有默认值; 静态初始化会回填
        ``gen_no``、``s_load``、``nbus``。单位沿用原 MATLAB 实现,功角为 rad,
        功率与内电势为标幺值。
    参数:
        m:机组惯性系数向量。
        d:机组阻尼系数向量。
        pmpu:机械功率(pu)。
        xd1:暂态电抗(pu)。
        epu:发电机内电势幅值(pu)。
        其余为故障线路、故障位置、ZIP 负荷比例等实验开关。
    """

    m: np.ndarray
    d: np.ndarray
    pmpu: np.ndarray
    xd1: np.ndarray
    epu: np.ndarray
    path_energy_cal: int = 0
    equ_cal: int = 2
    flag_xd: int = 0
    fault_line: np.ndarray = field(default_factory=lambda: np.array([9, 6], dtype=int))
    fault_position: int = 0
    p_load_zip: np.ndarray = field(default_factory=lambda: np.array([1.0, 0.0, 0.0]))
    q_load_zip: np.ndarray = field(default_factory=lambda: np.array([1.0, 0.0, 0.0]))
    gen_no: Optional[np.ndarray] = None
    s_load: Optional[np.ndarray] = None
    nbus: Optional[int] = None

    @property
    def ngen(self) -> int:
        # 用法:取发电机台数(惯量向量 m 的长度)。
        return int(self.m.size)

    @property
    def flag_uniform(self) -> bool:
        # 用法:判断是否为均匀阻尼(各机 d/m 相等),供能量法选择路径积分方式。
        ratios = self.d / self.m
        return bool(np.allclose(ratios, ratios[0], rtol=0.0, atol=1e-12))


@dataclass
class NetworkState:
    """单一网络工况的导纳、约简矩阵和 SEP 容器。

    使用方法:
        由网络构造流程填入 ``yfull`` 及 Kron 分块矩阵,求 SEP 后回填
        ``sep_delta``、``sep_omegapu``、``sep_perr``; ``metadata`` 存放重排、
        原始 yorg 等中间量供结构保持模型使用。
    参数说明:
        yred 为约简导纳; ynn/ynr/yrn/yrr 为发电机与负荷节点分块; 
        removed_bus 为被删除的故障母线编号(若有)。
    """

    name: str
    yfull: np.ndarray
    yred: np.ndarray
    ynn: np.ndarray
    ynr: np.ndarray
    yrn: np.ndarray
    yrr: np.ndarray
    removed_bus: Optional[int] = None
    sep_delta: Optional[np.ndarray] = None
    sep_omegapu: Optional[float] = None
    sep_perr: Optional[np.ndarray] = None
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class Trajectory:
    """约简模型固定步长 RK4 轨迹容器。

    使用方法:
        由 ``integrate_reduced`` 返回; ``theta``/``omega`` 为绝对量,
        ``thetac``/``omegac`` 为 COI 相对量,``pe`` 为逐步电磁功率,
        ``tunit`` 为步长(s)。
    """

    time: np.ndarray
    theta: np.ndarray
    omega: np.ndarray
    thetac: np.ndarray
    omegac: np.ndarray
    omegacoi: np.ndarray
    pe: np.ndarray
    tunit: float


@dataclass
class StaticResult:
    """对应 ``Cal_MM_Static`` 的一次静态初始化打包结果。

    使用方法:
        由 ``build_static_result`` 返回,作为后续实验(能量法、CCT、SPM)的统一
        输入; 内含 preset、basevalue、pfdata 及 pre/fault/post 三个网络工况。
    """

    preset: Preset
    basevalue: BaseValue
    pfdata: PFData
    prefault: NetworkState
    fault: NetworkState
    postfault: NetworkState
    emf: np.ndarray
    case: CaseData


@dataclass
class CCTResult:
    """能量法或时域搜索得到的临界切除时间结果。

    使用方法:
        读取 ``cct``(秒)与 ``flag_cct``(是否检测到越界); ``exit_state`` 存放
        对应时刻的 thetac/omegac,``method`` 标注来源方法,``detail`` 存放能量等
        中间量供检查。
    """

    cct: float
    exit_index: int
    exit_state: Dict[str, np.ndarray]
    flag_cct: bool
    method: str
    detail: Dict[str, Any] = field(default_factory=dict)
