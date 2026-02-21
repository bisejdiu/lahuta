# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         def __init__(self, v): self.v = v
#         def __or__(self, f): return Email(f(self.v))
#         def __str__(self): return self.v
#     print(str(
#         Email("")
#         | (lambda s: s + "besian")
#         | (lambda s: s + "sejdiu")
#         | (lambda s: s + "@gmail.com")
#     ))
#
"""Tests Conformer position operations using property-based testing."""

import numpy as np
import pytest

from lahuta import rdkit


# fmt: off
def test_hypothesis_conformer_positions():
    hyp = pytest.importorskip("hypothesis",             reason="Hypothesis not installed")
    st  = pytest.importorskip("hypothesis.strategies",  reason="Hypothesis not installed")
    hnp = pytest.importorskip("hypothesis.extra.numpy", reason="Hypothesis not installed")

    floats = st.floats(allow_infinity=False, allow_nan=False, width=64, min_value=-1e6, max_value=1e6)
    vec3 = hnp.arrays(np.float64, 3, elements=floats)
    n_strategy = st.integers(min_value=1, max_value=8)
    mat_n3 = n_strategy.flatmap(lambda n: hnp.arrays(np.float64, (n, 3), elements=floats))

    @hyp.given(vec3, vec3)
    def check_setAtomPos(v0: np.ndarray, v1: np.ndarray) -> None:
        c = rdkit.Conformer(2)
        c.setAtomPos(0, v0)
        c.setAtomPos(1, rdkit.Point3D(v1))
        p0 = c.getAtomPos(0)
        p1 = c.getAtomPos(1)
        assert np.allclose(np.array([p0.x, p0.y, p0.z]), v0)
        assert np.allclose(np.array([p1.x, p1.y, p1.z]), v1)

    @hyp.given(mat_n3)
    def check_setAllAtomPositions(mtx: np.ndarray) -> None:
        n = int(mtx.shape[0])
        c = rdkit.Conformer(n)
        c.setAllAtomPositions(mtx)
        got = c.getPositions()
        assert got.shape == (n, 3)
        assert np.allclose(got, mtx)

    check_setAtomPos()
    check_setAllAtomPositions()
