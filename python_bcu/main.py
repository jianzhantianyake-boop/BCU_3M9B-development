# -*- coding: utf-8 -*-
r"""BCU 3M9B Python 平台 - 中央控制台(一个入口跑通所有实验).

============================ 详细中文操作说明 ============================

[这是什么]
    整个 Python 平台的"总开关".菜单按 MATLAB 的 EXPERIMENT_MODE 组织成 8 个实验
    模式(8 个全部已实现), 再加几个工具项.
    你不用记 import, 运行本文件, 输入编号即可.

[运行前提]
    1) 已装 Python 3.12 与 numpy; 画图需要 matplotlib(已装).
    2) 在本文件所在目录 python_bcu\ 下运行(脚本会自动加好包路径).

[怎么运行] -- 在 PowerShell(提示符 PS C:\...> )里输入:

        cd python_bcu
        python main.py                 # 启动菜单, 按提示输入编号
        python main.py 3               # 直接跑模式 3(reduced_region)后退出
        python main.py a               # 依次跑模式 1-6
        python main.py help            # 打印本说明

[菜单 - 实验模式(对应 MATLAB EXPERIMENT_MODE)]
    1 reduced_cct        网络约简: 初始化 -> CCT/CUEP  (+ δ2-δ3 相平面图)
    2 reduced_numerical  网络约简: -> 三段数值轨迹      (2x2 多图)
    3 reduced_region     网络约简: 二维稳定域搜索        (平衡点 + 分界线)
    4 spm_cct            结构保持: 时域 CCT             (故障段约简/故障后 SPM)
    5 spm_numerical      结构保持: -> 数值轨迹          (3 图)
    6 spm_region         结构保持: 二维稳定域            (网格分类 + 解析分界线)
    7 two_machine_region_3d      两机完整模型: 三维稳定域 + 稳定流形曲面
    8 two_machine_region_3d_gfl  同上, 改用 GFL 参数(阻抗含电阻, 低惯量)

[菜单 - 工具项]
    p 交流潮流     s 静态初始化     t 冒烟自检
    a 依次跑模式 1-6     0 退出

[图片保存在哪]
    所有模式的图都存到 python_bcu\figures\ 目录(文件名即模式名).

[想改实验参数]
    下面 CONFIG 字典集中放了稳定域网格密度等常用旋钮; 更细的参数(时长/步长)在
    bcu_3m9b\experiments.py 各 mode_* 函数签名里, 改默认值即可.

[常见疑问]
    - "跑了没反应": 赋值/导入不打印; 本控制台每项都写了 print, 会显示结果与图路径.
    - "from ... 报错": 那是把 Python 代码敲进了 PowerShell; 本文件用 python main.py 运行.
    - 两机三维/GFL(模式 7/8): 已实现; 它们独立于 9 母线电网, 用两机完整模型 f_2m,
      GFL 只是换一组参数(阻抗含电阻, 低惯量), 不是另一套动力学.

=======================================================================
"""

from pathlib import Path
import sys

# 步骤0-A: Windows 控制台切到 UTF-8, 避免中文乱码.
if sys.platform == "win32":
    try:
        import ctypes
        ctypes.windll.kernel32.SetConsoleOutputCP(65001)
        ctypes.windll.kernel32.SetConsoleCP(65001)
    except Exception:
        pass
    for _stream in (sys.stdout, sys.stderr):
        try:
            _stream.reconfigure(encoding="utf-8")
        except Exception:
            pass

# 步骤0-B: 把本目录(包根)加入 sys.path, 任何工作目录下都能 import bcu_3m9b.
ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np


# ------------------------- 可调参数(改这里即可) -------------------------
# 用法: 稳定域网格越大越精细也越慢.
CONFIG = {
    "region_grid": 21,       # reduced_region 网格密度(每维点数)
    "spm_region_grid": 15,   # spm_region 网格密度(每维点数, 较慢)
}


# ------------------------- 实验模式(调用 experiments 模块) -------------------------

# 编号 -> (模式名, 说明).
MODE_ITEMS = [
    ("1", "reduced_cct", "网络约简: 初始化 -> CCT/CUEP (+ δ2-δ3 相平面图)"),
    ("2", "reduced_numerical", "网络约简: -> 三段数值轨迹 (2x2 多图)"),
    ("3", "reduced_region", "网络约简: 二维稳定域 (平衡点 + 分界线)"),
    ("4", "spm_cct", "结构保持: 时域 CCT (故障段约简/故障后 SPM)"),
    ("5", "spm_numerical", "结构保持: -> 数值轨迹 (3 图)"),
    ("6", "spm_region", "结构保持: 二维稳定域 (网格分类 + 解析分界线)"),
    ("7", "two_machine_region_3d", "两机完整模型: 三维稳定域 + 稳定流形曲面"),
    ("8", "two_machine_region_3d_gfl", "同上, GFL 参数 (阻抗含电阻, 低惯量)"),
]
PLACEHOLDER_ITEMS = []
_MODE_BY_KEY = {k: name for k, name, _ in MODE_ITEMS}


def run_mode(key: str, static=None) -> None:
    """按编号运行一个实验模式; 稳定域模式从 CONFIG 取网格密度.

    用法: 传入 '1'..'6'; 图会存到 figures/, 结果打印到控制台.
    """

    from bcu_3m9b import experiments

    name = _MODE_BY_KEY[key]
    func = experiments.MODES[name]
    kwargs = {}
    if name == "reduced_region":
        kwargs["grid_points"] = CONFIG["region_grid"]
    if name == "spm_region":
        kwargs["grid_points"] = CONFIG["spm_region_grid"]
    print(f"\n===== [{key}] {name} =====")
    func(static, **kwargs)


def show_placeholder(key: str) -> None:
    """打印未实现模式的说明(两机三维 / GFL)."""

    name = dict((k, n) for k, n, _ in PLACEHOLDER_ITEMS)[key]
    print(f"\n===== [{key}] {name} =====")
    print("该模式本次按你的选择未实现(只补单机与 SPM 六个模式).")
    print("需要时可后续补: 两机三维稳定域(Statable_Region_3D)或 GFL 版(Statable_Region_3D_GFL).")


# ------------------------- 工具项 -------------------------

def util_powerflow() -> None:
    """[p] 交流潮流: 求 case9 与 case39 的潮流解并打印残差."""

    from bcu_3m9b import case9_v2, case39_modified, solve_power_flow, to_pfdata

    print("\n===== [p] 交流潮流 =====")
    for name, case in (("case9_v2", case9_v2()), ("case39_modified", case39_modified())):
        pf = solve_power_flow(case, tol=1e-8)
        pfd = to_pfdata(pf)
        vmag = pfd.voltage[:, 0]
        print(f"[{name}] 收敛={pf.success} 迭代={pf.iterations} "
              f"残差范数={pf.residual_norm:.3e} "
              f"电压幅值范围=[{vmag.min():.4f}, {vmag.max():.4f}] pu")


def util_static() -> None:
    """[s] 静态初始化: 潮流 -> 网络约简 -> SEP, 打印关键量."""

    from bcu_3m9b import build_static_result

    print("\n===== [s] 静态初始化 =====")
    s = build_static_result()
    print(f"案例: {s.case.name}")
    print(f"母线={s.pfdata.nbus} 发电机={s.pfdata.ngen} 负荷={s.pfdata.nload}")
    print("预故障 SEP 角度 (rad):", np.array2string(s.prefault.sep_delta, precision=6))
    print("故障后 SEP 角度 (rad):", np.array2string(s.postfault.sep_delta, precision=6))
    print("故障后 SEP 残差:", np.array2string(s.postfault.sep_perr, precision=2))


def util_smoke() -> None:
    """[t] 冒烟自检: 调用 tests/smoke_test, 通过打印 PASS."""

    print("\n===== [t] 冒烟自检 =====")
    sys.path.insert(0, str(ROOT / "tests"))
    import smoke_test
    smoke_test.main()


UTILS = {"p": util_powerflow, "s": util_static, "t": util_smoke}


# ------------------------- 菜单与调度 -------------------------

def print_menu() -> None:
    """打印中文主菜单(8 实验模式 + 工具项)."""

    print("\n" + "=" * 60)
    print("   BCU 3M9B Python 平台 - 中央控制台")
    print("=" * 60)
    print("  实验模式(对应 MATLAB EXPERIMENT_MODE):")
    for key, name, desc in MODE_ITEMS:
        print(f"   [{key}] {name:<28} {desc}")
    for key, name, desc in PLACEHOLDER_ITEMS:
        print(f"   [{key}] {name:<28} {desc}")
    print("  工具项:")
    print("   [p] 交流潮流   [s] 静态初始化   [t] 冒烟自检")
    print("   [a] 依次跑模式 1-6      [0] 退出")
    print("=" * 60)


def dispatch(choice: str, static=None) -> None:
    """按输入分发到模式/占位/工具; 出错只打印简短原因."""

    try:
        if choice in _MODE_BY_KEY:
            run_mode(choice, static)
        elif choice in dict((k, n) for k, n, _ in PLACEHOLDER_ITEMS):
            show_placeholder(choice)
        elif choice in UTILS:
            UTILS[choice]()
        else:
            print(f"无效输入: {choice!r}")
    except Exception as exc:  # noqa: BLE001
        print(f"\n[运行 {choice} 时出错] {type(exc).__name__}: {exc}")
        print("提示: 确认在 python_bcu\\ 目录下运行, 且已装 numpy/matplotlib.")


def run_all() -> None:
    """依次跑模式 1-6.

    注意: 每个模式各自新做一次静态初始化, 不共用 static. 因为 run_bcu_experiment
    在求 CUEP 时会改写 postfault 的 SEP 字段, 共用会污染后续模式.
    """

    for key, _, _ in MODE_ITEMS:
        dispatch(key, None)


def interactive_loop() -> None:
    """交互式菜单主循环: 反复显示菜单, 读取编号, 运行, 直到输入 0."""

    while True:
        print_menu()
        try:
            choice = input("请输入编号后回车(0 退出): ").strip().lower()
        except (EOFError, KeyboardInterrupt):
            print("\n已退出.")
            return
        if choice == "0":
            print("已退出.")
            return
        if choice == "a":
            run_all()
        else:
            dispatch(choice)


def main(argv=None) -> None:
    """入口: 无参数进菜单; 带参数直接运行.

    用法:
        python main.py            进入交互菜单
        python main.py 3          直接跑模式 3 后退出
        python main.py p          跑工具项(交流潮流)
        python main.py a          依次跑模式 1-6
        python main.py help       打印顶部说明
    """

    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        interactive_loop()
        return
    arg = argv[0].strip().lower()
    if arg in ("help", "-h", "--help"):
        print(__doc__)
    elif arg == "a":
        run_all()
    else:
        dispatch(arg)


if __name__ == "__main__":
    main()
