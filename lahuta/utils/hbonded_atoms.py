"""
Module: hbonded_atoms.py

This module provides a utility function to find the hydrogen-bonded atoms in a molecule. It leverages 
MDAnalysis and OpenBabel libraries to accomplish this.

Functions:
    find_hydrogen_bonded_atoms(mda, mol): Finds the hydrogen-bonded 
                                            atoms in the molecule.

Notes:
    The function 'find_hydrogen_bonded_atoms' first establishes a mapping between the atoms in the given 
    AtomGroup and their respective indices. Then, it iterates over each atom in the molecule and checks if it 
    has any bonded hydrogen atoms. If yes, it updates the corresponding row in the hydrogen bond array with 
    the indices of the bonded hydrogen atoms.
"""

import numpy as np

# from numpy.typing import NDArray
from openbabel import openbabel as ob
from scipy.sparse import csr_matrix, dok_matrix

from lahuta.lahuta_types.openbabel import MolType


def find_hydrogen_bonded_atoms(mol: MolType, n_atoms: int) -> csr_matrix:
    """
    Identifies the hydrogen-bonded atoms in a molecule.

    This function first establishes a mapping between the atoms in the given AtomGroup and their respective indices.
    Then, it iterates over each atom in the molecule and checks if it has any bonded hydrogen atoms. If yes, it
    updates the corresponding row in the hydrogen bond array with the indices of the bonded hydrogen atoms.

    Args:
        mda (AtomGroupType): An AtomGroup from MDAnalysis, representing the atoms in a molecule.
        mol (MolType): An OpenBabel OBMol object, representing the molecule.

    Returns:
        NDArray[np.int32]: A 2D numpy array of shape (n_atoms, 6) where each row corresponds to an atom in the molecule
                           and contains the indices of the hydrogen atoms bonded to it. A value of -1 indicates no
                           bonded hydrogen atom for that position.
    """
    hbond_array = dok_matrix((n_atoms, 6), dtype=np.int32)

    for atom in ob.OBMolAtomIter(mol):
        if atom.ExplicitHydrogenCount():
            for ix, atom2 in enumerate(ob.OBAtomAtomIter(atom)):
                if atom2.GetAtomicNum() == 1:
                    atom1_id = atom.GetId()
                    atom2_id = atom2.GetId()
                    hbond_array[atom1_id, ix] = atom2_id

    return hbond_array.tocsr()
