"""脚本用途：运行一个完整的 BCU/能量/CCT/暂态实验。

使用方法：
    命令行执行 ``python scripts/run_experiment.py``，无需参数；脚本先做静态初始化
    再跑非交互式实验并打印判据。
"""

from pathlib import Path
import sys

# 步骤1：把包根目录加入 sys.path，免安装即可导入。
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from bcu_3m9b.bcu import build_static_result, run_bcu_experiment


def main() -> None:
    """执行非交互式实验并打印判据。

    步骤：
        (1) 构造静态结果；(2) 以给定故障时长/步长跑实验；(3) 打印故障轨迹、退出点、
        MGP、CUEP 来源、临界能量和 LEA/REA CCT。
    """

    static = build_static_result()
    result = run_bcu_experiment(
        static, fault_time=0.2, tunit=1e-3, postfault_time=2.0, cct_samples=11
    )
    # 打印实验结果（输出全英文）。
    print("=== BCU Python experiment platform ===")
    print(f"Fault trajectory length: {result['fault_trajectory'].time[-1]:.6g} s")
    print(f"Exit point: index={result['exit_index']}, t={result['exit_time']:.6g} s")
    print("MGP:", result["mgp"]["theta_mgp"])
    print("CUEP source:", result["cuep_source"])
    print("Critical energy:", f"{result['critical_energy']:.12g}")
    print("LEA CCT:", f"{result['lea'].cct:.6g} s, flag={result['lea'].flag_cct}")
    print("REA CCT grid estimate:", f"{result['rea_cct']:.6g} s")
    print("Note: CUEP/MGP and CCT are currently runnable experiment interfaces; strict cross-validation is established separately.")


if __name__ == "__main__":
    main()
