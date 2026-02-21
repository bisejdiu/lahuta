# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         a = "besian"
#         b = "sejdiu"
#         c = "@gmail.com"
#     print(getattr(Email, "a") + getattr(Email, "b") + getattr(Email, "c"))
#
"""Tests atom queries, PDB residue info, and pooled atom info."""

import pytest

from lahuta import rdkit


def test_atom_queries_basic():
    m = rdkit.RWMol()
    i = m.addAtom()
    a = m.getAtomWithIdx(i)
    a.setAtomicNum(6)  # Carbon

    # For a default "empty" atom, check basic invariants
    assert a.getIdx() == i
    assert a.getDegree() == 0
    assert a.getFormalCharge() == 0
    assert a.getIsAromatic() in (True, False)

    a.calcImplicitValence()
    assert a.getNumExplicitHs() >= 0 and a.getNumImplicitHs() >= 0
    assert isinstance(a.getMass(), float)
    assert isinstance(a.getSymbol(), str) and len(a.getSymbol()) >= 1

    hyb = a.getHybridization()
    assert hasattr(hyb, "name") and hasattr(hyb, "value")
    assert a.getMonomerInfo() is None  # No monomer info by default


def test_atom_pdb_residue_info_construction_and_setters():
    info = rdkit.AtomPDBResidueInfo(
        atomName=" CA ",
        serialNumber=10,
        residueName="GLY",
        residueNumber=7,
        chainId="A",
    )
    assert info.getSerialNumber() == 10
    assert info.getResidueName() == "GLY"
    assert info.getResidueNumber() == 7
    assert info.getChainId() == "A"

    # defs
    assert info.getAltLoc() == ""
    assert info.getInsertionCode() == ""
    assert info.getOccupancy() == pytest.approx(1.0)
    assert info.getTempFactor() == pytest.approx(0.0)
    assert info.getSegmentNumber() == 0
    assert info.getIsHeteroAtom() is False
    assert info.getSecondaryStructure() == 0

    # setters
    info.setAltLoc("B")
    info.setInsertionCode(" ")
    info.setOccupancy(0.75)
    info.setTempFactor(12.5)
    info.setIsHeteroAtom(True)
    info.setSecondaryStructure(2)
    info.setSegmentNumber(3)
    info.setResidueIndex(5)
    info.setResidueName("ALA")
    info.setResidueNumber(8)
    info.setChainId("Z")

    assert info.getAltLoc() == "B"
    assert info.getInsertionCode() == " "
    assert info.getOccupancy() == pytest.approx(0.75)
    assert info.getTempFactor() == pytest.approx(12.5)
    assert info.getSegmentNumber() == 3
    assert info.getResidueIndex() == 5
    assert info.getResidueName() == "ALA"
    assert info.getResidueNumber() == 8
    assert info.getChainId() == "Z"
    assert info.getIsHeteroAtom() is True
    assert info.getSecondaryStructure() == 2


def test_pooled_atom_pdb_residue_info_overloads_and_reset():
    pinfo = rdkit.pAtomPDBResidueInfo()
    # Initialize later
    pinfo.initialize(" N  ", 1, "MET", 1)
    assert pinfo.getResidueName() == "MET" and pinfo.getResidueNumber() == 1

    # Construct via overload with common fields
    pinfo2 = rdkit.pAtomPDBResidueInfo(" O  ", 2, "MET", 1)
    assert pinfo2.getSerialNumber() == 2
    assert pinfo2.getResidueName() == "MET"

    # mutate then reset
    pinfo2.setAltLoc("")
    pinfo2.setChainId("A")
    pinfo2.setIsHeteroAtom(True)
    pinfo2.resetState()

    assert pinfo2.getAltLoc() == ""
    assert pinfo2.getChainId() == "A"  # NOTE: we don't reset chainId
    # rdkit has this False by default, but our pooled extension doesn't reset it
    assert pinfo2.getIsHeteroAtom() is True
