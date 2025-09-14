"""Provides a utility function to find the hydrogen-bonded atoms in a molecule."""

import numpy as np
from openbabel import openbabel as ob
from scipy.sparse import csr_matrix, dok_matrix

from lahuta._types.openbabel import MolType


def find_hydrogen_bonded_atoms(mol: MolType, n_atoms: int) -> csr_matrix:
    """Identify the hydrogen-bonded atoms in a molecule.

    This function takes a molecule and the number of atoms in that molecule as input, and it
    returns a sparse matrix with the indices of the hydrogen-bonded atoms.

    Args:
        mol (MolType): The molecule for which the hydrogen-bonded atoms are to be identified.
        n_atoms (int): The number of atoms in the molecule.

    Returns:
        (csr_matrix): A sparse matrix with the indices of the hydrogen-bonded atoms.
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
