"""
Module: radii.py

This module contains utility functions for determining the van der Waals of atoms in a molecule.
It provides two functions that achieve the same goal but in different ways; the first one 
iterates over the atoms in the molecule which can be time-consuming, 
while the second one utilizes a vectorized operation for faster performance.

Functions:
    assign_radii(mol): Determines the radii for each atom in the molecule using iteration.

    v_radii_assignment(elements): Determines the radii for each atom 
                                    in the molecule using vectorized operation.

"""
from typing import Dict

import numpy as np
from MDAnalysis.topology.tables import vdwradii as MDA_VDW_RADII  # type: ignore
from numpy.typing import NDArray
from openbabel import openbabel as ob

from lahuta.lahuta_types.openbabel import MolType


def assign_radii(mol: MolType) -> NDArray[np.float32]:
    """Get the van der Waals for each atom in the molecule.

    Parameters
    ----------
    mol : openbabel.OBMol
        The molecule object.

    Returns
    -------
    radii_array:  An array of shape (n_atoms,) containing the van der Waals
                    radii for each atom in the molecule.
    """
    atom_radii = np.zeros(mol.NumAtoms(), dtype=float)
    for atom in ob.OBMolAtomIter(mol):
        atom_radii[atom.GetId()] = ob.GetVdwRad(atom.GetAtomicNum())  # type: ignore
    return atom_radii


def v_radii_assignment(
    elements: NDArray[np.str_],
) -> NDArray[np.float32]:
    """
    Determines the van der Waals radii for each atom in the molecule using vectorized operation.

    This function is a faster alternative to 'assign_radii' as it utilizes numpy's vectorized operations
    for improved performance when dealing with large molecules.

    Args:
        elements: An array of atomic element symbols.

    Returns:
        NDArray[np.float32]: An array of shape (n_atoms, ) where each element is the van der Waals radius.
    """

    vdwradii: Dict[str, int] = {k.capitalize(): v for k, v in MDA_VDW_RADII.items()}  # type: ignore

    def v_capitalize(array: NDArray[np.str_], mapping: Dict[str, int]) -> NDArray[np.float32]:
        vfunc = np.vectorize(mapping.get)
        return vfunc(array)  # type: ignore

    return v_capitalize(elements, vdwradii)
