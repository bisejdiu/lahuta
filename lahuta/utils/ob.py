"""
Utility functions for OpenBabel.
"""

from typing import List

import numpy as np
from numpy.typing import NDArray
from openbabel import openbabel as ob

from lahuta.lahuta_types.openbabel import BondIterable, MolType, ObRingType, ObVector3Wrapper


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


class Rings:
    def __init__(self) -> None:
        self._centers: List[NDArray[np.float32]] = []
        self._normals: List[NDArray[np.float32]] = []
        self._atoms: List[List[int]] = []
        self._first_atom_idx: List[int] = []

    def add_ring(self, ob_ring: ObRingType) -> None:
        ob_vector3_wrapper = ObVector3Wrapper(ob.vector3(), ob.vector3())
        center = ob_vector3_wrapper.center
        normal = ob_vector3_wrapper.normal
        ob_ring.findCenterAndNormal(center, normal, ob.vector3())
        center_coords: NDArray[np.float32] = np.array([center.GetX(), center.GetY(), center.GetZ()])
        normal_coords: NDArray[np.float32] = np.array([normal.GetX(), normal.GetY(), normal.GetZ()])
        atoms = sorted([atom for atom in ob_ring._path])  # type: ignore

        self._centers.append(center_coords)
        self._normals.append(normal_coords)
        self._atoms.append(atoms)
        self._first_atom_idx.append(atoms[0])

    @property
    def centers(self) -> NDArray[np.float32]:
        return np.array(self._centers, dtype=np.float32)

    @property
    def normals(self) -> NDArray[np.float32]:
        return np.array(self._normals, dtype=np.float32)

    @property
    def atoms(self) -> NDArray[np.int32]:
        return np.array(self._atoms, dtype=object)

    @property
    def first_atom_idx(self) -> NDArray[np.int32]:
        return np.array(self._first_atom_idx, dtype=np.int32)

    def __len__(self) -> int:
        return len(self._centers)


def enumerate_rings(mol: MolType) -> Rings:
    rings = Rings()
    for ob_ring in mol.GetSSSR():
        if ob_ring.IsAromatic():
            rings.add_ring(ob_ring)
    return rings
