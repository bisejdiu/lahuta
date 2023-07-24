"""
Placeholder for the neighbors module.
"""

from typing import Any, Dict, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from lahuta.config.defaults import CONTACTS, VDW_RADII
from lahuta.core.helpers import get_class_attributes
from lahuta.lahuta_types.mdanalysis import AtomGroupType
from lahuta.lahuta_types.openbabel import MolType
from lahuta.utils import array_utils as au
from lahuta.utils.atom_types import find_hydrogen_bonded_atoms
from lahuta.utils.math import calc_pairwise_distances, calc_vertex_angles
from lahuta.writers.frame_writer import DataFrameWriter


class HBondHandler:
    def __init__(self, atoms: AtomGroupType, hbond_array: NDArray[np.int32]):
        self._atoms = atoms
        self.hbond_array = hbond_array

    def get_hbond_distances(self, attr_col: AtomGroupType, hbound_attr_col: AtomGroupType) -> NDArray[np.float32]:
        hbound_atom_indices = self.hbond_array[hbound_attr_col.atoms.indices]
        hbound_atom_pos = self._atoms.positions[hbound_atom_indices]

        hbound_atom_pos[hbound_atom_indices == 0] = np.nan
        distance_array = calc_pairwise_distances(attr_col.atoms.positions, hbound_atom_pos)

        return distance_array

    def get_vdw_distances(self, attr_col: AtomGroupType, vdw_comp_factor: float) -> NDArray[np.float32]:
        return attr_col.atoms.vdw_radii + VDW_RADII["H"] + vdw_comp_factor

    def get_hbond_angles(self, col1: AtomGroupType, col2: AtomGroupType) -> NDArray[np.float32]:
        atom1_pos = col1.atoms.positions
        atom2_pos = col2.atoms.positions

        hbound_atom_indices = self.hbond_array[col1.atoms.indices]
        hbound_atom_pos = self._atoms.positions[hbound_atom_indices]

        hbound_atom_pos[hbound_atom_indices == 0] = np.nan

        point_a = atom1_pos[:, np.newaxis, :]
        point_b = hbound_atom_pos
        point_c = atom2_pos[:, np.newaxis, :]

        return calc_vertex_angles(point_b, point_a, point_c, degrees=False)


class NeighborPairs:
    """A class for storing neighbor pairs."""

    def __init__(
        self,
        mda: AtomGroupType,
        mol: MolType,
        pairs: NDArray[np.int32],
        distances: NDArray[np.float32],
    ):
        self.mda = mda
        self.mol = mol
        self.atoms = self.mda.atoms.universe.atoms

        self._validate_inputs(pairs, distances)
        self._pairs, self._distances = NeighborPairs.sort_inputs(pairs, distances)

        self.hbond_array: NDArray[np.int32] = find_hydrogen_bonded_atoms(self.mda, self.mol)
        self.hbond_handler = HBondHandler(self.atoms, self.hbond_array)
        self.hbond_angles: NDArray[np.float32] = np.array([])
        self._annotations: Dict[str, NDArray[Any]] = {}

    def _validate_inputs(self, pairs: NDArray[np.int32], distances: NDArray[np.float32]) -> None:
        message = (
            "The number of pairs and distances must be the same."
            f"Got {pairs.shape[0]} pairs and {distances.shape[0]} distances."
        )
        assert pairs.shape[0] == distances.shape[0], message

    @staticmethod
    def get_sorting_index(pairs: NDArray[np.int32]) -> NDArray[np.int32]:
        sorted_pairs: NDArray[np.int32] = np.sort(pairs, axis=1)
        indices: NDArray[np.int32] = np.argsort(sorted_pairs[:, 0])  # type: ignore

        return indices

    @staticmethod
    def sort_inputs(
        pairs: NDArray[np.int32], distances: NDArray[np.float32]
    ) -> Tuple[NDArray[np.int32], NDArray[np.float32]]:
        pairs = np.sort(pairs, axis=1)
        indices = np.argsort(pairs[:, 0])

        return pairs[indices], distances[indices]

    def _get_pair_column(self, partner: int) -> AtomGroupType:
        """Return the column of the pair of atoms depending on the value of partner."""
        return self.atoms[self.pairs[:, partner - 1]]

    def _get_partners(self, partner: int) -> Tuple[AtomGroupType, AtomGroupType]:
        """Return attr_col and hbound_attr_col depending on the value of partner."""
        partner2 = 1 if partner == 2 else 2

        return self._get_pair_column(partner), self._get_pair_column(partner2)

    def type_filter(self, atom_type: str, partner: int) -> "NeighborPairs":
        """Filters pairs based on atom types.

        Args:
            atom_type (str or list of str): The atom types can one of the following:
                'carbonyl_oxygen', 'weak_hbond_donor', 'pos_ionisable', 'carbonyl_carbon',
                'hbond_acceptor', 'hbond_donor', 'neg_ionisable', 'weak_hbond_acceptor',
                'xbond_acceptor', 'aromatic', 'hydrophobe'.
                Names come from :class:`SmartsPatternRegistry` (Enum).

            partner (int): The column to select the atom types from. It can be either 1 or 2.

        Returns:
            NeighborPairs: A NeighborPairs object containing the filtered pairs.
        """
        col_ag = getattr(self._get_pair_column(partner), atom_type)
        mask = col_ag.astype(bool)

        return self.clone(self.pairs[mask], self.distances[mask])

    def index_filter(
        self,
        indices: NDArray[np.int32],
        partner: int,
    ) -> "NeighborPairs":
        """Select pairs based on the atom indices.

        Parameters
        ----------
        indices : int or list of int
            The atom indices to select.

        partner : int
            The column to select the atom types from. Either 1 or 2.

        Returns
        -------
        pairs : NeighborPairs
            A NeighborPairs object containing the selected pairs.
        """

        col_func = self._get_pair_column(partner)
        mask = np.isin(col_func.indices, indices)  # type: ignore

        return self.clone(self.pairs[mask], self.distances[mask])

    def distance_filter(self, distance: float) -> "NeighborPairs":
        """Select pairs based on the distance.

        Parameters
        ----------
        distance : float
            The distance to select.

        Returns
        -------
        pairs : NeighborPairs
            A NeighborPairs object containing the selected pairs.
        """
        mask = self.distances <= distance
        return self.clone(self.pairs[mask], self.distances[mask])

    def numeric_filter(self, array: NDArray[np.float32], cutoff: float, lte: bool = True) -> "NeighborPairs":
        """Select pairs based on a boolean mask.

        Parameters
        ----------
        mask : np.ndarray
            The boolean mask to select.

        Returns
        -------
        pairs : NeighborPairs
            A NeighborPairs object containing the selected pairs.
        """
        # add support for lt and gt
        mask = array <= cutoff if lte else array > cutoff
        return self.clone(self.pairs[mask], self.distances[mask])

    def radius_filter(self, radius: float, partner: int) -> "NeighborPairs":
        """Select pairs based on the radius.

        Parameters
        ----------
        radius : float
            The radius to select.

        partner : int
            The column to select the atom types from. Either 1 or 2.

        Returns
        -------
        pairs : NeighborPairs
            A NeighborPairs object containing the selected pairs.
        """

        col_func = self._get_pair_column(partner)
        mask = col_func.atoms.vdw_radii <= radius

        return self.clone(self.pairs[mask], self.distances[mask])

    def hbond_distance_filter(self, partner: int, vdw_comp_factor: float = 0.1) -> "NeighborPairs":
        """Filter the pairs based on the distance between the hydrogen bonded atoms.

        Parameters
        ----------
        partner : int
            The column of the hydrogen bonded atom indices in the `hbond_array`.
        vdw_comp_factor : float
            The van der Waals complementarity factor.

        Returns
        -------
        pairs : NeighborPairs
            The filtered neighbor pairs.
        """
        attr_col, hbound_attr_col = self._get_partners(partner)

        vdw_distances = self.hbond_handler.get_vdw_distances(attr_col, vdw_comp_factor)
        hbond_dist = self.hbond_handler.get_hbond_distances(attr_col, hbound_attr_col)

        distances_mask = np.any(hbond_dist <= vdw_distances[:, np.newaxis], axis=1)  # type: ignore
        hbond_dist_pairs = self.pairs[distances_mask]
        hbond_distances = self.distances[distances_mask]

        return self.clone(hbond_dist_pairs, hbond_distances)

    def hbond_angle_filter(self, partner: int, weak: bool = False) -> "NeighborPairs":
        """Filter the pairs based on the angle between the hydrogen bonded atoms.

        Parameters
        ----------
        partner : int
            The column of the hydrogen bonded atom indices in the `hbond_array`.
        weak : bool
            If True, accept weaker hydrogen bonds.

        Returns
        -------
        pairs : NeighborPairs
            The filtered neighbor pairs.
        """
        contact_type = "weak hbond" if weak else "hbond"
        attr_partner, hbound_attr_partner = self._get_partners(partner)

        # if self.hbond_angles is None:
        self.hbond_angles = self.hbond_handler.get_hbond_angles(attr_partner, hbound_attr_partner)

        idx = np.any(self.hbond_angles >= CONTACTS[contact_type]["angle rad"], axis=1)  # type: ignore
        self._pairs = self._pairs[idx]
        self._distances = self._distances[idx]

        return self.clone(self.pairs, self.distances)

    def intersection(self, other: "NeighborPairs") -> "NeighborPairs":
        """Return the intersection of two NeighborPairs objects.

        Parameters
        ----------
        other : NeighborPairs
            The other NeighborPairs object.

        Returns
        -------
        pairs : NeighborPairs
            A NeighborPairs object containing the intersection of the two NeighborPairs objects.
        """

        mask = au.intersection(self.pairs, other.pairs)

        return self.clone(self.pairs[mask], self.distances[mask])

    def union(self, other: "NeighborPairs") -> "NeighborPairs":
        """Return the union of two NeighborPairs objects.

        Parameters
        ----------
        other : NeighborPairs
            The other NeighborPairs object.

        Returns
        -------
        pairs : NeighborPairs
            A NeighborPairs object containing the union of the two NeighborPairs objects.
        """
        mask = au.union(self.pairs, other.pairs)
        pairs = np.concatenate((self.pairs, other.pairs), axis=0)[mask]  # type: ignore
        distances = np.concatenate((self.distances, other.distances), axis=0)[mask]  # type: ignore

        return self.clone(pairs, distances)

    def difference(self, other: "NeighborPairs") -> "NeighborPairs":
        """Return the difference of two NeighborPairs objects.

        Parameters
        ----------
        other : NeighborPairs
            The other NeighborPairs object.

        Returns
        -------
        pairs : NeighborPairs
            A NeighborPairs object containing the difference of the two NeighborPairs objects.
        """
        mask = au.difference(self.pairs, other.pairs)

        return self.clone(self.pairs[mask], self.distances[mask])

    def symmetric_difference(self, other: "NeighborPairs") -> "NeighborPairs":
        """Return the symmetric difference of two NeighborPairs objects.

        Parameters
        ----------
        other : NeighborPairs
            The other NeighborPairs object.

        Returns
        -------
        pairs : NeighborPairs
            A NeighborPairs object containing the symmetric difference of the two NeighborPairs objects.
        """
        mask_a, mask_b = au.symmetric_difference(self.pairs, other.pairs)

        pairs = np.concatenate((self.pairs[mask_a], other.pairs[mask_b]), axis=0)  # type: ignore
        distances = np.concatenate((self.distances[mask_a], other.distances[mask_b]), axis=0)  # type: ignore

        unique_indices = np.unique(pairs, axis=0, return_index=True)[1]  # type: ignore
        sorted_indices = np.sort(unique_indices)  # type: ignore

        pairs = pairs[sorted_indices]
        distances = distances[sorted_indices]

        sorted_pairs = pairs[np.argsort(pairs[:, 0])]  # type: ignore
        sorted_distances = distances[np.argsort(pairs[:, 0])]  # type: ignore

        return self.clone(sorted_pairs, sorted_distances)

    def isdisjoint(self, other: "NeighborPairs") -> bool:
        """Test whether two NeighborPairs objects have a null intersection.

        Parameters
        ----------
        other : NeighborPairs
            The other NeighborPairs object.

        Returns
        -------
        isdisjoint : bool
            True if the two NeighborPairs objects have a null intersection.
        """
        return au.isdisjoint(self.pairs, other.pairs)

    def issubset(self, other: "NeighborPairs") -> bool:
        """Test whether all elements of a NeighborPairs object are in another.

        Parameters
        ----------
        other : NeighborPairs
            The other NeighborPairs object.

        Returns
        -------
        issubset : bool
            True if all elements of the NeighborPairs object are in the other NeighborPairs object.
        """
        return au.issubset(self.pairs, other.pairs)

    def issuperset(self, other: "NeighborPairs") -> bool:
        """Test whether all elements of another NeighborPairs object are in this one.

        Parameters
        ----------
        other : NeighborPairs
            The other NeighborPairs object.

        Returns
        -------
        issuperset : bool
            True if all elements of the other NeighborPairs object are in this NeighborPairs object.
        """
        return au.issuperset(self.pairs, other.pairs)

    def isequal(self, other: "NeighborPairs") -> bool:
        """Test whether two NeighborPairs objects contain the same elements.

        Parameters
        ----------
        other : NeighborPairs
            The other NeighborPairs object.

        Returns
        -------
        isequal : bool
            True if the two NeighborPairs objects contain the same elements.
        """
        return au.isequal(self.pairs, other.pairs)

    def isunique(self) -> bool:
        """Test whether the NeighborPairs object contains unique pairs.

        Returns
        -------
        isunique : bool
            True if the NeighborPairs object contains unique pairs.
        """
        return au.isunique(self.pairs)

    def is_strict_subset(self, other: "NeighborPairs") -> bool:
        """Test whether all elements of a NeighborPairs object are in another and the two sets are not equal.

        Parameters
        ----------
        other : NeighborPairs
            The other NeighborPairs object.

        Returns
        -------
        is_strict_subset : bool
            True if all elements of the NeighborPairs object are in the other NeighborPairs object and the two sets are not equal.
        """
        return au.is_strict_subset(self.pairs, other.pairs)

    def is_strict_superset(self, other: "NeighborPairs") -> bool:
        """Test whether all elements of another NeighborPairs object are in this one and the two sets are not equal.

        Parameters
        ----------
        other : NeighborPairs
            The other NeighborPairs object.

        Returns
        -------
        is_strict_superset : bool
            True if all elements of the other NeighborPairs object are in this NeighborPairs object and the two sets are not equal.
        """
        return au.is_strict_superset(self.pairs, other.pairs)

    def clone(self, pairs: NDArray[np.int32], distances: NDArray[np.float32]) -> "NeighborPairs":
        """Get a copy of the NeighborPairs object."""

        attrs = {attr: getattr(self, attr) for attr in get_class_attributes(self)}
        attrs.update({"_pairs": pairs, "_distances": distances})

        cls = type(self)
        child_instance = cls.__new__(cls)

        for attr, value in attrs.items():
            if not isinstance(getattr(cls, attr, None), property):
                setattr(child_instance, attr, value)

        return child_instance

    @property
    def annotations(self) -> Dict[str, NDArray[Any]]:
        """Get the annotations of the NeighborPairs object."""
        return self._annotations

    @annotations.setter
    def annotations(self, annotations: Dict[str, NDArray[Any]]) -> None:
        """Set the annotations of the NeighborPairs object."""
        self._annotations = annotations

    def add_annotations(self, annotations: Dict[str, NDArray[Any]]) -> None:
        """Add annotations to the NeighborPairs object."""
        for value in annotations.values():
            assert len(value) == self.pairs.shape[0]

        self._annotations.update(annotations)

    def to_frame(
        self,
        df_format: Literal["compact", "expanded"] = "expanded",
        annotations: bool = False,
    ) -> pd.DataFrame:
        """Convert the NeighborPairs object to a pandas DataFrame.

        Parameters
        ----------
        df_format : str
            The format of the DataFrame. It can be either "compact" or "expanded".
        annotations : bool
            Whether to include annotations in the DataFrame.

        Returns
        -------
        df : pd.DataFrame
            The DataFrame containing the pairs of atoms and their distances.
        """
        if annotations:
            return self._create_df(df_format, self.annotations)
        else:
            return self._create_df(df_format)

        # return self._create_df(df_format)

    def to_dict(self, df_format: Literal["compact", "expanded"] = "expanded") -> Dict[str, Any]:
        """Convert the NeighborPairs object to a dictionary.

        Parameters
        ----------
        df_format : str
            The format of the DataFrame. It can be either "compact" or "expanded".

        Returns
        -------
        df : pd.DataFrame
            The DataFrame containing the pairs of atoms and their distances.
        """
        return self._create_df(df_format).to_dict(orient="list")  # type: ignore

    def _create_df(
        self,
        df_format: Literal["compact", "expanded"] = "expanded",
        annotations: Optional[Dict[str, NDArray[Any]]] = None,
    ) -> pd.DataFrame:
        """Create a DataFrame from the NeighborPairs object."""
        return DataFrameWriter(self, df_format, annotations).create()

    def _neighborpairs_equal(self, other: "NeighborPairs") -> bool:
        # Get the indices that would sort each array
        indices1: NDArray[np.int32] = np.lexsort((self.pairs[:, 1], self.pairs[:, 0]))  # type: ignore
        indices2: NDArray[np.int32] = np.lexsort((other.pairs[:, 1], other.pairs[:, 0]))  # type: ignore

        # Sort each array using the indices
        pairs, dists = self.pairs[indices1], self.distances[indices1]
        other_pairs, other_dists = other.pairs[indices2], other.distances[indices2]

        return np.array_equal(pairs, other_pairs) and np.array_equal(dists, other_dists)  # type: ignore

    @property
    def partner1(self) -> AtomGroupType:
        """Get the first partner of the pairs of atoms that are neighbors."""
        return self._get_pair_column(1)

    @property
    def partner2(self) -> AtomGroupType:
        """Get the second partner of the pairs of atoms that are neighbors."""
        return self._get_pair_column(2)

    @property
    def pairs(self) -> NDArray[np.int32]:
        """Get the pairs of atoms that are neighbors."""
        return self._pairs

    @property
    def distances(self) -> NDArray[np.float32]:
        """Get the distances between the pairs of atoms that are neighbors."""
        return self._distances

    @property
    def indices(self) -> NDArray[np.int32]:
        """Get the indices of the atoms that are neighbors."""
        return np.unique([self.partner1.indices, self.partner2.indices])  # type: ignore

    def __getitem__(self, item: Union[int, slice, NDArray[np.int32]]) -> "NeighborPairs":
        """Get the pair of atoms at the specified index.

        Parameters
        ----------
        item : int or slice
            The index of the pair of atoms to get.

        Returns
        -------
        NeighborPair
            A new NeighborPair object containing the pair of atoms at the specified index.
        """
        if isinstance(item, int):
            return self.clone(
                self.pairs[item],
                self.distances[item],
            )

        # return self.__class__(self._atoms, self.pairs[item], self.distances[item])
        return self.clone(self.pairs[item], self.distances[item])

    def __contains__(self, other: "NeighborPairs") -> bool:
        """Test whether a pair of atoms is in the NeighborPairs object.

        Parameters
        ----------
        other : NeighborPair
            The pair of atoms to test.

        Returns
        -------
        contains : bool
            True if the pair of atoms is in the NeighborPairs object.
        """

        return au.issubset(other.pairs, self.pairs)

    def __add__(self, other: "NeighborPairs") -> "NeighborPairs":
        # TODO:
        # currently this is different from MDAnalysis.
        # The question to answer is if neighbor pairs should be unique or not.
        return self.union(other)

    def __sub__(self, other: "NeighborPairs") -> "NeighborPairs":
        return self.difference(other)

    def __or__(self, other: "NeighborPairs") -> "NeighborPairs":
        return self.symmetric_difference(other)

    def __eq__(self, other: Any) -> bool:
        return self._neighborpairs_equal(other)

    def __and__(self, other: "NeighborPairs") -> "NeighborPairs":
        return self.intersection(other)

    def __xor__(self, other: "NeighborPairs") -> "NeighborPairs":
        return self.symmetric_difference(other)

    def __len__(self) -> int:
        return self.pairs.shape[0]

    def __str__(self) -> str:
        return f"<Lahuta NeighborPairs class containing {self.indices.size} atoms and {self.pairs.shape[0]} pairs>"

    def __repr__(self) -> str:
        return self.__str__()
