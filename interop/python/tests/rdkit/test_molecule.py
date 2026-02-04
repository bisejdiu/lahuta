# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     def gen():
#         yield from ["besian", "sejdiu", "@gmail.com"]
#     print("".join(gen()))
#
"""Tests molecule construction, atoms, bonds, and queries."""

from lahuta import rdkit


def test_rw_mol_atoms_bonds_and_queries():
    m = rdkit.RWMol()

    a0 = m.addAtom()
    a1 = m.addAtom()
    a2 = m.addAtom()
    assert (a0, a1, a2) == (0, 1, 2)

    assert m.getNumAtoms() == 3
    assert m.getNumAtoms(True) == 3
    assert m.getNumBonds() == 0
    assert m.getNumBonds(True) == 0

    nb = m.addBond(0, 1, rdkit.BondType.SINGLE)
    assert nb == 1
    nb = m.addBond(1, 2, rdkit.BondType.DOUBLE)
    assert nb == 2
    assert m.getNumBonds() == 2

    # bond pairs and objects
    pairs = m.bondPairs()
    assert set(map(tuple, pairs)) == {(0, 1), (1, 2)}
    b_list = m.bonds()
    assert len(b_list) == 2
    assert all(isinstance(b, rdkit.Bond) for b in b_list)

    # neighbors and incident bonds
    assert set(m.atomNeighbors(1)) == {0, 2}
    assert set(map(tuple, m.atomBonds(1))) == {(0, 1), (1, 2)}
    assert len(m.atomBondObjects(1)) == 2

    # getBondBetweenAtoms and getBondWithIdx
    b01 = m.getBondBetweenAtoms(0, 1)
    assert isinstance(b01, rdkit.Bond)
    assert b01.getBeginAtomIdx() == 0 and b01.getEndAtomIdx() == 1
    assert b01.getOtherAtomIdx(0) == 1 and b01.getOtherAtomIdx(1) == 0

    b0 = m.getBondWithIdx(0)
    assert b0.getBondType() == rdkit.BondType.SINGLE
    # mutate via reference and observe on refetch
    b0.setIsAromatic(True)
    assert m.getBondWithIdx(0).getIsAromatic() is True

    # These don't work yet. At some point I'll fix them if necessary.
    # Set and get direction/stereo
    # b0.setBondDir(rd.BondDir.BEGINWEDGE)
    # assert b0.getBondDir() == rd.BondDir.BEGINWEDGE
    # b0.setStereo(rd.BondStereo.STEREOCIS)
    # assert b0.getStereo() == rd.BondStereo.STEREOCIS
    # # Stereo reference atoms (accepts 0 or 2 indices)
    # b0.setStereoAtoms(0, 2)
    # refs = b0.getStereoAtoms()
    # assert isinstance(refs, list) and len(refs) in (0, 2)

    # Atoms retrieval and iteration
    a = m.getAtomWithIdx(1)
    assert a.getIdx() == 1
    assert m[1].getIdx() == 1
    atoms_seq = list(m.__iter__())
    assert len(atoms_seq) == m.getNumAtoms()

    # Connected components: add two isolated atoms as a new components
    m.addAtom(6)
    m.addAtom(6)
    n_comp, mapping = m.getConnectedComponents()
    assert isinstance(n_comp, int) and isinstance(mapping, list)
    assert len(mapping) == m.getNumAtoms()
    # chain + two isolates = 3 components
    assert n_comp == 3

    # # Utilities: property caches and clear props
    m.updatePropertyCache()
    m.clearComputedProps()

    # Remove bond and atom
    m.removeBond(1, 2)
    assert m.getNumBonds() == 1
    m.removeAtom(2)
    assert m.getNumAtoms() == 4
