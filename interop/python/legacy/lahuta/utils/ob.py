"""Utility functions for OpenBabel."""

import numpy as np
from numpy.typing import NDArray
from openbabel import openbabel as ob

from lahuta._types.openbabel import MolType, ObRingType, ObVector3Wrapper


def get_bonded_atoms(mol: MolType) -> NDArray[np.int32]:
    """Return an array of bonded atoms in a molecule.

    Args:
        mol (MolType): The molecule for which to find bonded atoms.

    Returns:
        NDArray[np.int32]: A numpy array containing the indices of bonded atoms.
    """
    bonds: NDArray[np.int32] = np.zeros((mol.NumBonds(), 2), dtype=int)
    for ix, bond in enumerate(ob.OBMolBondIter(mol)):
        atom_idx1, atom_idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonds[ix, :] = (atom_idx1, atom_idx2) if atom_idx1 < atom_idx2 else (atom_idx2, atom_idx1)

    return bonds


class Rings:
    """A class for storing information about aromatic rings in a molecule.

    This class stores the center and normal vectors for each aromatic ring in a molecule.
    It also stores the indices of the atoms in each ring and the index of the first atom
    in each ring.
    """

    def __init__(self) -> None:
        self.rings = []
        self._centers: list[NDArray[np.float32]] = []
        self._normals: list[NDArray[np.float32]] = []
        self._atoms: list[list[int]] = []
        self._first_atom_idx: list[int] = []

    def add_ring(self, ob_ring: ObRingType) -> None:
        """Add a ring to the Rings object.

        Args:
            ob_ring (ObRingType): The ring to be added.
        """
        self.rings.append(ob_ring)
        ob_vector3_wrapper = ObVector3Wrapper(ob.vector3(), ob.vector3())
        center = ob_vector3_wrapper.center
        normal = ob_vector3_wrapper.normal
        ob_ring.findCenterAndNormal(center, normal, ob.vector3())
        center_coords: NDArray[np.float32] = np.array([center.GetX(), center.GetY(), center.GetZ()])
        normal_coords: NDArray[np.float32] = np.array([normal.GetX(), normal.GetY(), normal.GetZ()])
        # FIX: _path is not really the atom id?
        atoms = sorted([x - 1 for x in ob_ring._path])  # noqa: SLF001
        # atoms = sorted([x for x in ob_ring._path])  # noqa: SLF001

        self._centers.append(center_coords)
        self._normals.append(normal_coords)
        self._atoms.append(atoms)
        atoms.sort()
        self._first_atom_idx.append(atoms[0])

    @property
    def centers(self) -> NDArray[np.float32]:
        """Returns the centers of the rings.

        Returns:
            NDArray[np.float32]: The centers of the rings.
        """
        return np.array(self._centers, dtype=np.float32)

    @property
    def normals(self) -> NDArray[np.float32]:
        """Returns the normal vectors of the rings.

        Returns:
            NDArray[np.float32]: The normal vectors of the rings.
        """
        return np.array(self._normals, dtype=np.float32)

    @property
    def atoms(self) -> NDArray[np.int32]:
        """Returns the atoms in the rings.

        Returns:
            NDArray[np.int32]: The atoms in the rings.
        """
        return np.array(self._atoms, dtype=object)

    @property
    def first_atom_idx(self) -> NDArray[np.int32]:
        """Returns the index of the first atom in each ring.

        Returns:
            NDArray[np.int32]: The index of the first atom in each ring.
        """
        return np.array(self._first_atom_idx, dtype=np.int32)

    def __len__(self) -> int:
        return len(self._centers)


def enumerate_rings(mol: MolType) -> Rings:
    """Enumerate the aromatic rings in a molecule.

    This function takes a molecule as input and returns a Rings object containing information
    about the aromatic rings in the molecule. It is used by `atom-plane` and `plane-plane`
    contact calculations.

    Args:
        mol (MolType): The molecule for which to enumerate the rings.

    Returns:
        Rings: A Rings object containing information about the rings in the molecule.
    """
    rings = Rings()
    for ob_ring in mol.GetSSSR():
        if ob_ring.IsAromatic():
            rings.add_ring(ob_ring)
    return rings
