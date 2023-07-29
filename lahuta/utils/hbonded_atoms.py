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
    # n_atoms: int = mda.universe.atoms.n_atoms
    # print('NATOMS', n_atoms)
    # hbond_array: NDArray[np.int32] = np.zeros((n_atoms, 6), dtype=int)
    hbond_array = dok_matrix((n_atoms, 6), dtype=np.int32)

    # will give -1 for atoms not in the atomgroup
    # TODO: & FIXME: We need to use or re-use csc or dok matrix
    # max_index = np.max(mda.indices)
    # atom_mapping = np.full(max_index + 1, -1)
    # atom_mapping[np.arange(mda.n_atoms)] = mda.indices
    # print("ATOM MAPPING DETAILS: ", atom_mapping.shape)

    for atom in ob.OBMolAtomIter(mol):
        if atom.ExplicitHydrogenCount():
            for ix, atom2 in enumerate(ob.OBAtomAtomIter(atom)):  # type: ignore
                if atom2.GetAtomicNum() == 1:  # type: ignore
                    atom1_id = atom.GetId()  # type: ignore
                    atom2_id = atom2.GetId()  # type: ignore
                    hbond_array[atom1_id, ix] = atom2_id

    return hbond_array.tocsr()
