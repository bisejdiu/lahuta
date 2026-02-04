# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     print(functools.reduce(operator.add, ["besian", "sejdiu", "@gmail.com"]))
#
"""Tests Conformer construction, position management, and basic operations."""

import numpy as np
import pytest

from lahuta import rdkit


def test_conformer_basic_and_positions_numpy_and_point3d():

    c = rdkit.Conformer(2)
    assert c.getNumAtoms() == 2
    assert isinstance(c.getId(), int)
    c.setId(42)
    assert c.getId() == 42
    assert c.is3D() in (True, False)
    c.set3D(True)
    assert c.is3D() is True

    # Positions view
    pos = c.getPositions()
    assert isinstance(pos, np.ndarray)
    assert pos.shape == (2, 3)
    # Default-initialized to zeros
    assert np.allclose(pos, 0.0)

    # setAtomPos via numpy
    c.setAtomPos(0, np.array([1.0, 2.0, 3.0], dtype=np.float64))
    p0 = c.getAtomPos(0)
    assert (pytest.approx(p0.x) == 1.0) and (pytest.approx(p0.y) == 2.0) and (pytest.approx(p0.z) == 3.0)

    # setAtomPos via Point3D
    p1 = rdkit.Point3D(-1.0, -2.0, -3.0)
    c.setAtomPos(1, p1)
    q1 = c.getAtomPos(1)
    assert (q1.x, q1.y, q1.z) == (-1.0, -2.0, -3.0)

    # setAllAtomPositions
    all_pos = np.array([[0.0, 0.0, 0.5], [4.0, 5.0, 6.0]], dtype=np.float64)
    c.setAllAtomPositions(all_pos)
    assert np.allclose(c.getPositions(), all_pos)

    # Resize and reserve should behave
    old_n = c.getNumAtoms()
    c.reserve(10)
    assert c.getNumAtoms() == old_n  # reserve doesn't change size
    c.resize(3)
    assert c.getNumAtoms() == 3
    assert c.getPositions().shape == (3, 3)


def test_conformer_view_is_mutable_and_reflects_back():
    c = rdkit.Conformer(1)
    arr = c.getPositions()
    assert arr.flags["C_CONTIGUOUS"]
    arr[0, :] = [7.0, 8.0, 9.0]
    p = c.getAtomPos(0)
    assert (p.x, p.y, p.z) == (7.0, 8.0, 9.0)


def test_has_non_zero_z_coords():
    c = rdkit.Conformer(3)
    c.set3D(True)
    # all z == 0 -> False
    assert rdkit.hasNonZeroZCoords(c) is False
    # set a non-zero z
    c.setAtomPos(1, np.array([0.0, 0.0, 1.0], dtype=np.float64))
    assert rdkit.hasNonZeroZCoords(c) is True
