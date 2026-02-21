# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     a, *rest = ["besian", "sejdiu", "@gmail.com"]
#     print(a + "".join(rest))
#
import numpy as np
import pytest

from lahuta import rdkit as lrdkit


def test_conformer_positions_zero_copy_view_and_mutation_roundtrip():
    conf = lrdkit.Conformer(3)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float64)

    conf.setAllAtomPositions(coords)
    arr = conf.getPositions()

    assert arr.shape == (3, 3)
    assert arr.dtype == np.float64
    assert arr.flags["OWNDATA"] is False
    assert arr.base is conf

    # Mutate numpy view, Conformer should reflect the change
    x0 = conf.getAtomPos(0).x
    arr[0, 0] = x0 + 0.5
    assert conf.getAtomPos(0).x == x0 + 0.5

    # Mutate Conformer, numpy view should reflect the change
    conf.setAtomPos(1, np.array([9.0, 8.0, 7.0], dtype=np.float64))
    assert np.allclose(arr[1], np.array([9.0, 8.0, 7.0], dtype=np.float64))


def test_conformer_set_positions_shape_errors():
    conf = lrdkit.Conformer(2)
    with pytest.raises(ValueError):
        conf.setAllAtomPositions(np.ones((2, 2), dtype=np.float64))  # wrong columns
    with pytest.raises(ValueError):
        conf.setAtomPos(0, np.ones((2,), dtype=np.float64))  # shape (3,) required
