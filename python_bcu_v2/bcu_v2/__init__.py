# -*- coding: utf-8 -*-
"""bcu_v2: 在冻结的 v1(bcu_3m9b)之上做正确性验证与数值增强.

设计:
    v1(../python_bcu/bcu_3m9b)只读复用; 本包只新增 P0(验证套件 + 隐患修复)与
    P1.2(事件驱动的精确 CCT). 导入本包时会把兄弟目录 python_bcu 加入 sys.path.
"""

from __future__ import annotations

import sys
from pathlib import Path

# 路径引导: 把兄弟目录 python_bcu 加入 sys.path, 使 `import bcu_3m9b` 可用.
_V1 = Path(__file__).resolve().parents[2] / "python_bcu"
if _V1.is_dir() and str(_V1) not in sys.path:
    sys.path.insert(0, str(_V1))

__all__ = ["smib", "cct", "fixes", "invariants"]
