"""A module containing the HBondHandler class.

Classes:
    HBondHandler: A class used to compute various properties of hydrogen bonds for a given atomic group.

"""
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    from lahuta import NeighborPairs
import numpy as np
from numpy.typing import NDArray
from scipy.sparse import coo_array, csr_matrix

from lahuta._types.mdanalysis import AtomGroupType
from lahuta.config.defaults import CONTACTS, VDW_RADII
from lahuta.utils.hbonded_atoms import find_hydrogen_bonded_atoms
from lahuta.utils.math import calc_pairwise_distances, calc_vertex_angles


class HBondHandler:
    """Compute various properties of hydrogen bonds for a given atomic group.

    Attributes
    ----------
    _atoms : AtomGroupType
        The group of atoms under consideration.
    hbond_array : NDArray[np.int32]
        The indices of the hydrogen bonded atoms.
    """

    def __init__(self, ns: "NeighborPairs"):
        self.atoms = ns.atoms
        self.hbond_array = find_hydrogen_bonded_atoms(ns.luni.to("mol"), self.atoms.n_atoms)

    def get_hbond_distances(self, attr_col: AtomGroupType, hbound_attr_col: AtomGroupType) -> NDArray[np.float32]:
        """Compute the distances between hydrogen atoms and their respective bonded atoms.

        Args:;
            attr_col (AtomGroupType):
                Group of atoms for which the distances are computed.
            hbound_attr_col (AtomGroupType):
                Group of atoms to which the hydrogen atoms are bonded.

        Returns:
            NDArray[np.float32]:
                An array of the distances between each atom in attr_col and its corresponding
                bonded atom in hbound_attr_col.
        """
        indices = hbound_attr_col.atoms.indices
        hbond_array = cast(csr_matrix, self.hbond_array[indices])
        selected_rows_coo: coo_array = hbond_array.tocoo()
        hbound_atom_indices = np.zeros_like(selected_rows_coo.toarray(), dtype=np.int32)
        hbound_atom_indices[selected_rows_coo.row, selected_rows_coo.col] = selected_rows_coo.data

        hbound_atom_pos = np.take(self.atoms.positions, hbound_atom_indices, axis=0)
        hbound_atom_pos[hbound_atom_indices == 0] = np.nan

        return calc_pairwise_distances(attr_col.atoms.positions, hbound_atom_pos)

    def get_vdw_distances(self, attr_col: AtomGroupType, vdw_comp_factor: float) -> NDArray[np.float32]:
        """Compute the Van der Waals distances between the atoms in the given group.

        Args:
            attr_col (AtomGroupType):
                Group of atoms for which the Van der Waals distances are computed.
            vdw_comp_factor (float):
                A factor to adjust the computation of Van der Waals distances.

        Returns:
            NDArray[np.float32]:
                An array of the Van der Waals distances for each atom in attr_col.
        """
        return attr_col.atoms.vdw_radii + VDW_RADII["H"] + vdw_comp_factor

    def get_hbond_angles(self, col1: AtomGroupType, col2: AtomGroupType) -> NDArray[np.float32]:
        """Compute the angles formed by the hydrogen bonds between atoms in col1 and col2.

        Args:
            col1 (AtomGroupType):
                Group of atoms that form one end of the bond.
            col2 (AtomGroupType):
                Group of atoms that form the other end of the bond.

        Returns:
            NDArray[np.float32]:
                An array of the angles (in radians) for each hydrogen bond between atoms in col1 and col2.
        """
        atom1_pos = col1.atoms.positions
        atom2_pos = col2.atoms.positions

        indices = col1.atoms.indices
        hbond_array = cast(csr_matrix, self.hbond_array[indices])
        selected_rows_coo: coo_array = hbond_array.tocoo()
        hbound_atom_indices = np.zeros_like(selected_rows_coo.toarray(), dtype=np.int32)
        hbound_atom_indices[selected_rows_coo.row, selected_rows_coo.col] = selected_rows_coo.data

        hbound_atom_pos = np.take(self.atoms.positions, hbound_atom_indices, axis=0)
        hbound_atom_pos[hbound_atom_indices == 0] = np.nan

        point_a = atom1_pos[:, np.newaxis, :]
        point_b = hbound_atom_pos
        point_c = atom2_pos[:, np.newaxis, :]

        return calc_vertex_angles(point_b, point_a, point_c, degrees=False)
    
    def hbond_distance_filter(self, ns: "NeighborPairs", partner: int, vdw_comp_factor: float = 0.1) -> "NeighborPairs":
        """Filter the pairs based on the distance between the hydrogen bonded atoms.

        The method filters pairs from the NeighborPairs object where the hydrogen bond distances
        are less than or equal to the specified van der Waals distances. The `partner` parameter specifies
        the column of hydrogen bonded atom indices in the `hbond_array`.

        Args:
            ns (NeighborPairs): A NeighborPairs object containing the atom neighbor relationships in the system.
            partner (int): The column of the hydrogen bonded atom indices in the `hbond_array`.
            vdw_comp_factor (float, optional): The van der Waals complementarity factor. Defaults to 0.1.

        Returns:
            A NeighborPairs object containing the pairs that meet the hydrogen bond distance filter.
        """
        attr_col, hbound_attr_col = ns.get_partners(partner)

        vdw_distances = self.get_vdw_distances(attr_col, vdw_comp_factor)
        hbond_dist = self.get_hbond_distances(attr_col, hbound_attr_col)

        distances_mask = np.any(hbond_dist <= vdw_distances[:, np.newaxis], axis=1)
        hbond_dist_pairs = ns.pairs[distances_mask]
        hbond_distances = ns.distances[distances_mask]

        return ns.new(hbond_dist_pairs, hbond_distances)

    def hbond_angle_filter(self, ns: "NeighborPairs", partner: int, weak: bool = False) -> "NeighborPairs":
        """Filter the pairs based on the angle between the hydrogen bonded atoms.

        The method filters pairs from the NeighborPairs object where the hydrogen bond angles are greater
        than or equal to the specified contact angle. The `partner` parameter specifies the column of hydrogen
        bonded atom indices in the `hbond_array`. If `weak` is True, the function will accept weaker hydrogen bonds.

        Args:
            ns (NeighborPairs): A NeighborPairs object containing the atom neighbor relationships in the system.
            partner (int): The column of the hydrogen bonded atom indices in the `hbond_array`.
            weak (bool, optional): If True, accept weaker hydrogen bonds. Defaults to False.

        Returns:
            A NeighborPairs object containing the pairs that meet the hydrogen bond angle filter.
        """
        contact_type = "weak hbond" if weak else "hbond"
        attr_partner, hbound_attr_partner = ns.get_partners(partner)

        self.hbond_angles = self.get_hbond_angles(attr_partner, hbound_attr_partner)

        idx = np.any(self.hbond_angles >= CONTACTS[contact_type]["angle rad"], axis=1)
        pairs = ns.pairs[idx]
        distances = ns.distances[idx]

        return ns.new(pairs, distances)
