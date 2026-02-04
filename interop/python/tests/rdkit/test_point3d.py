# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     T = type("", (), {"__getattr__": lambda self, name: {
#         "f": "besian", "l": "sejdiu", "d": "@gmail.com"
#     }[name]})
#     print(T().f + T().l + T().d)
#
"""Tests Point3D construction, string representation, and numpy array integration."""

import numpy as np

from lahuta import rdkit


def test_point3d_construction_and_str():
    p0 = rdkit.Point3D()
    assert p0.x == 0.0 and p0.y == 0.0 and p0.z == 0.0
    s = str(p0)
    assert s.startswith("(") and s.endswith(")") and "," in s

    p1 = rdkit.Point3D(1.0, -2.5, 3.25)
    assert (p1.x, p1.y, p1.z) == (1.0, -2.5, 3.25)

    arr = np.array([9.0, 8.0, 7.0], dtype=np.float64)
    p2 = rdkit.Point3D(arr)
    assert (p2.x, p2.y, p2.z) == (9.0, 8.0, 7.0)
