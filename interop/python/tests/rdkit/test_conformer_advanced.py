"""Tests Conformer ownership, copying, and molecule integration."""

import numpy as np


def test_conformer_mol_ownership_and_copying(rdkit):
    m = rdkit.RWMol()
    c = rdkit.Conformer(2)
    # Add matching atoms to the molecule so conformer can be attached
    for _ in range(c.getNumAtoms()):
        m.addAtom()
    # addConformer copies, original should not be owned by a mol
    assert c.hasOwningMol() is False
    cid = m.addConformer(c, assign_id=True)
    assert isinstance(cid, int)
    assert m.getNumConformers() == 1
    mc = m.getConformer(cid)
    assert mc.hasOwningMol() is True
    assert np.allclose(mc.getPositions(), c.getPositions())
    # removal works
    m.removeConformer(cid)
    assert m.getNumConformers() == 0
