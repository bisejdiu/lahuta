# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     print("".__class__("besian") + "".__class__("sejdiu") + "".__class__("@gmail.com"))
#
"""Tests that enums are accessible and work correctly with bond operations."""

from lahuta import rdkit


def test_enums_accessible_and_roundtrip_on_bond():
    m = rdkit.RWMol()
    m.addAtom()
    m.addAtom()
    m.addBond(0, 1, rdkit.BondType.TRIPLE)
    b = m.getBondWithIdx(0)
    assert b.getBondType() == rdkit.BondType.TRIPLE

    for d in (rdkit.BondDir.NONE, rdkit.BondDir.BEGINWEDGE, rdkit.BondDir.BEGINDASH):
        b.setBondDir(d)
        assert b.getBondDir() == d
    for s in (rdkit.BondStereo.STEREONONE, rdkit.BondStereo.STEREOANY):
        b.setStereo(s)
        assert b.getStereo() == s
