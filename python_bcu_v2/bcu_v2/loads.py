# -*- coding: utf-8 -*-
"""P2.3: ZIP 负荷模型(恒阻抗 Z / 恒电流 I / 恒功率 P 组合).

使用方法:
    zip_load_power(V, P0, Q0, zp, zq) 给出电压相关的负荷功率; zip_algebraic_residual 把 ZIP
    负荷接入结构保持模型的负荷母线功率平衡(推广 v1 只支持恒功率的代数残差).

约定:
    zp=[aZ, aI, aP], 三者和为 1: P(V)=P0*(aZ*V^2 + aI*V + aP). Q 同理用 zq.
    极限: aP=1 恒功率(= v1 行为); aZ=1 恒阻抗; aI=1 恒电流.

说明: 本模块提供 ZIP 的模型函数与残差, 可插入 spm_dae 的代数求解; 完整多母线联立留待接线.
"""

from __future__ import annotations

import numpy as np


def zip_load_power(V, P0, Q0, zp=(0.0, 0.0, 1.0), zq=(0.0, 0.0, 1.0)):
    """计算 ZIP 负荷在电压 V 下的有功/无功.

    使用方法: 传入电压幅值 V(标量或数组)、额定 P0/Q0、ZIP 系数 zp/zq(各三元, 和为1);
    返回 (P, Q). 默认恒功率.
    """

    V = np.asarray(V, dtype=float)
    zp = np.asarray(zp, dtype=float)
    zq = np.asarray(zq, dtype=float)
    fp = zp[0] * V ** 2 + zp[1] * V + zp[2]
    fq = zq[0] * V ** 2 + zq[1] * V + zq[2]
    return P0 * fp, Q0 * fq


def zip_algebraic_residual(z, delta_gen, yorg_ordered, load_pq0, ngen,
                           zp=(0.0, 0.0, 1.0), zq=(0.0, 0.0, 1.0)):
    """ZIP 负荷下的结构保持负荷母线功率平衡残差(推广 v1.algebraic_residual).

    使用方法: z=[负荷角(nload); 负荷电压(nload)]; load_pq0 为额定 [P0, Q0]; 返回残差向量.
    ZIP 令负荷功率随母线电压变化: 网络注入 + P_zip(V) = 0, 无功同理.
    """

    z = np.asarray(z, dtype=float)
    nload = yorg_ordered.shape[0] - ngen
    delta_load = z[:nload]
    voltage_load = z[nload:]
    voltage = np.r_[np.exp(1j * np.asarray(delta_gen)),
                    voltage_load * np.exp(1j * delta_load)]
    injection = voltage * np.conj(yorg_ordered @ voltage)
    P_zip, Q_zip = zip_load_power(voltage_load, load_pq0[:, 0], load_pq0[:, 1], zp, zq)
    return np.r_[injection.real[ngen:] + P_zip, injection.imag[ngen:] + Q_zip]
