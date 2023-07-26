"""
Utility functions for OpenBabel.
"""

import numpy as np
from numpy.typing import NDArray
from openbabel import openbabel as ob

from lahuta.lahuta_types.openbabel import BondIterable, MolType


def get_bonded_atoms(mol: MolType) -> NDArray[np.int32]:
    """
    Returns an array of bonded atoms in a molecule.

    Args:
        mol (MolType): The molecule for which to find bonded atoms.

    Returns:
        NDArray[np.int32]: A numpy array containing the indices of bonded atoms.
    """

    def bond_iter_wrapper(mol: MolType) -> BondIterable:
        """Wrapper for the openbabel bond iterator."""
        return ob.OBMolBondIter(mol)  # type: ignore

    bonds: NDArray[np.int32] = np.zeros((mol.NumBonds(), 2), dtype=int)
    for ix, bond in enumerate(bond_iter_wrapper(mol)):
        atom_idx1, atom_idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonds[ix, :] = (atom_idx1, atom_idx2) if atom_idx1 < atom_idx2 else (atom_idx2, atom_idx1)

    return bonds
