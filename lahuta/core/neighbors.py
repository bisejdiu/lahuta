"""
Placeholder for the neighbors module.
"""

import warnings
from abc import abstractmethod
from functools import partial, update_wrapper
from typing import Any, List, Optional, Union

import numpy as np
from typing_extensions import Literal, Protocol

from ..config.defaults import CONTACTS, VDW_RADII
from ..utils import array_utils as au
from ..utils.array_utils import calculate_angle
from ..utils.atom_types import find_hydrogen_bonded_atoms


class NeighborPairsBase(Protocol):
    """A protocol for neighbor pairs."""

    _atoms: Any
    col1: Any
    col2: Any
    _pairs: np.ndarray
    _distances: np.ndarray

    @abstractmethod
    def type_filter(
        self, atom_types: Union[str, List[str]], col: int
    ) -> "NeighborPairsBase":
        """Filter the neighbor pairs by atom type."""

    @abstractmethod
    def contact_type(
        self, contact_type: Literal["hbond", "native"]
    ) -> Union["HBondNeighborPairs", "NeighborPairs"]:
        """Filter the neighbor pairs by contact type."""

    @abstractmethod
    def index_filter(
        self, indices: Union[int, List[int]], col: int
    ) -> "NeighborPairsBase":
        """Filter the neighbor pairs by atom index."""

    @abstractmethod
    def distance_filter(self, distance: float) -> "NeighborPairsBase":
        """Filter the neighbor pairs by distance."""

    @abstractmethod
    def radius_filter(self, radius: float, col: int) -> "NeighborPairsBase":
        """Filter the neighbor pairs by radius."""

    @abstractmethod
    def hbond_distance_filter(
        self, col: int, vdw_comp_factor: float
    ) -> "HBondNeighborPairs":
        """Filter the neighbor pairs by distance."""

    @property
    def pairs(self) -> np.ndarray:
        """Return the neighbor pairs."""
        return self._pairs

    @property
    def distances(self) -> np.ndarray:
        """Return the distances of the neighbor pairs."""
        return self._distances

    @property
    def indices(self) -> np.ndarray:
        """Get the indices of the atoms that are neighbors."""
        return np.unique([self.col1.indices, self.col2.indices])


class NeighborPairs:
    """A class for storing neighbor pairs."""

    def __init__(self, uniatom, pairs, distances, **kwargs):
        self._atoms = uniatom.atoms
        self._pairs = pairs
        self._distances = distances
        self._angles = None
        self.col1, self.col2 = uniatom[pairs[:, 0]], uniatom[pairs[:, 1]]

        self.setops = au.ArraySetOps(self._pairs)

        # put kwargs in the object
        for key, value in kwargs.items():
            setattr(self, key, value)

        assert (
            self.col1.indices.size == self.col2.indices.size == self._distances.size
        ), "The first and second pairs columns must have the same size as the distances."

        # TODO: Get from config file and change variable name
        self.type_keys = {
            "hbond acceptor": 0,
            "pos ionisable": 1,
            "carbonyl oxygen": 2,
            "weak hbond donor": 3,
            "carbonyl carbon": 4,
            "weak hbond acceptor": 5,
            "hbond donor": 6,
            "neg ionisable": 7,
            "aromatic": 8,
            "xbond acceptor": 9,
            "hydrophobe": 10,
        }

    def contact_type(self, contact_type: Literal["hbond", "native"]):
        """Filter the `NeighborPairs` based on the contact type.

        Parameters
        ----------
        contact_type : str
            The contact type to filter by. Must be either 'hbond' or 'native'.

        Returns
        -------
        NeighborPairs or HBondNeighborPairs
            A `NeighborPairs` or `HBondNeighborPairs` object depending on the `contact_type` parameter.
        """
        if contact_type == "native":
            return self
        elif contact_type == "hbond":
            return HBondNeighborPairs(self._atoms, self._pairs, self._distances)
        else:
            raise ValueError(
                f"contact_type must be 'hbond' or 'native', not {contact_type}"
            )

    def type_filter(
        self, atom_types: Union[str, List[str]], col: int
    ) -> "NeighborPairs":
        """Select pairs based on the atom types.

        Parameters
        ----------
        atom_types : str or list of str
            The atom types to select. The atom types can be a combination of the
            following: 'carbonyl oxygen', 'weak hbond donor', 'pos ionisable',
            'carbonyl carbon', 'hbond acceptor', 'hbond donor', 'neg ionisable',
            'weak hbond acceptor', 'xbond acceptor', 'aromatic', 'hydrophobe'.

        col : AtomGroup
            The column to select the atom types from. Either 1 or 2.

        Returns
        -------
        pairs : NeighborPairs
            A NeighborPairs object containing the selected pairs.
        """
        if isinstance(atom_types, str):
            atom_types = [atom_types]

        col_func = getattr(self, f"col{col+1}")
        mask = np.any(
            col_func.atom_types[:, [self.type_keys[k] for k in atom_types]], axis=1
        )

        return self.__class__(self._atoms, self.pairs[mask], self.distances[mask])

    def index_filter(
        self,
        indices: Union[int, List[int]],
        col: int,
    ) -> Union["NeighborPairs", "HBondNeighborPairs"]:
        """Select pairs based on the atom indices.

        Parameters
        ----------
        indices : int or list of int
            The atom indices to select.

        col : int
            The column to select the atom types from. Either 1 or 2.

        Returns
        -------
        pairs : NeighborPairs
            A NeighborPairs object containing the selected pairs.
        """
        if isinstance(indices, int):
            indices = [indices]

        col_func = getattr(self, f"col{col+1}")
        mask = np.isin(col_func.indices, indices)

        return self.__class__(self._atoms, self.pairs[mask], self.distances[mask])

    def distance_filter(
        self, distance: float
    ) -> Union["NeighborPairs", "HBondNeighborPairs"]:
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
        return self.__class__(self._atoms, self.pairs[mask], self.distances[mask])

    def radius_filter(
        self, radius: float, col: int
    ) -> Union["NeighborPairs", "HBondNeighborPairs"]:
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

        col_func = getattr(self, f"col{col+1}")
        mask = col_func.atoms.vdw_radii <= radius

        return self.__class__(self._atoms, self.pairs[mask], self.distances[mask])

    def angle_filter(
        self, angle: float
    ) -> Union["NeighborPairs", "HBondNeighborPairs"]:
        """Select pairs based on the angle.

        Parameters
        ----------
        angle : float
            The angle to select.

        Returns
        -------
        pairs : NeighborPairs
            A NeighborPairs object containing the selected pairs.
        """
        if self._angles is None:

            warnings.warn("Found no angles. Returning all pairs.")
            return self

        mask = self._angles <= angle
        return self.__class__(self._atoms, self.pairs[mask], self.distances[mask])

    # TODO: vdw_comp_factor should be a class attribute and retrieved from config file
    def hbond_distance_filter(self, col: int = 0, vdw_comp_factor: float = 0.1):
        raise NotImplementedError(
            "This method is not implemented for `NeighborPairs`. Please use `HBondNeighborPairs`."
        )

    def intersection(
        self, other: "NeighborPairsBase"
    ) -> Union["NeighborPairs", "HBondNeighborPairs"]:
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

        mask = self.setops.intersection(other.pairs)

        return self.__class__(
            self._atoms,
            self.pairs[mask],
            self._distances[mask],
        )

    def union(
        self, other: "NeighborPairsBase"
    ) -> Union["NeighborPairs", "HBondNeighborPairs"]:
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
        mask = self.setops.union(other.pairs)
        pairs = np.concatenate((self.pairs, other.pairs), axis=0)[mask]
        distances = np.concatenate((self.distances, other.distances), axis=0)[mask]

        return self.__class__(
            self._atoms,
            pairs,
            distances,
        )

    def difference(
        self, other: "NeighborPairsBase"
    ) -> Union["NeighborPairs", "HBondNeighborPairs"]:
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
        mask = self.setops.difference(other.pairs)

        return self.__class__(
            self._atoms,
            self.pairs[mask],
            self._distances[mask],
        )

    def symmetric_difference(
        self, other: "NeighborPairsBase"
    ) -> Union["NeighborPairs", "HBondNeighborPairs"]:
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
        mask_a, mask_b = self.setops.symmetric_difference(other.pairs)

        pairs = np.concatenate((self.pairs[mask_a], other.pairs[mask_b]), axis=0)
        distances = np.concatenate(
            (self.distances[mask_a], other.distances[mask_b]), axis=0
        )

        unique_indices = np.unique(pairs, axis=0, return_index=True)[1]
        sorted_indices = np.sort(unique_indices)

        pairs = pairs[sorted_indices]
        distances = distances[sorted_indices]

        sorted_pairs = pairs[np.argsort(pairs[:, 0])]
        sorted_distances = distances[np.argsort(pairs[:, 0])]

        return self.__class__(
            self._atoms,
            sorted_pairs,
            sorted_distances,
        )

    def isdisjoint(self, other: "NeighborPairsBase") -> bool:
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
        return self.setops.isdisjoint(other.pairs)

    def issubset(self, other: "NeighborPairsBase") -> bool:
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
        return self.setops.issubset(other.pairs)

    def issuperset(self, other: "NeighborPairsBase") -> bool:
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
        return self.setops.issuperset(other.pairs)

    def isequal(self, other: "NeighborPairsBase") -> bool:
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
        return self.setops.isequal(other.pairs)

    def isunique(self) -> bool:
        """Test whether the NeighborPairs object contains unique pairs.

        Returns
        -------
        isunique : bool
            True if the NeighborPairs object contains unique pairs.
        """
        return self.setops.isunique()

    def is_strict_subset(self, other: "NeighborPairsBase") -> bool:
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
        return self.setops.is_strict_subset(other.pairs)

    def is_strict_superset(self, other: "NeighborPairsBase") -> bool:
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
        return self.setops.is_strict_superset(other.pairs)

    @property
    def pairs(self) -> np.ndarray:
        """Get the pairs of atoms that are neighbors."""
        return self._pairs

    @property
    def distances(self) -> np.ndarray:
        """Get the distances between the pairs of atoms that are neighbors."""
        return self._distances

    @property
    def angles(self) -> Optional[np.ndarray]:
        """Get the angles between the pairs of atoms that are neighbors."""
        return self._angles

    @property
    def indices(self) -> np.ndarray:
        """Get the indices of the atoms that are neighbors."""
        return np.unique([self.col1.indices, self.col2.indices])

    def __getitem__(self, item: Union[int, slice, np.ndarray]) -> "NeighborPairs":
        """Get the pair of atoms at the specified index.

        Parameters
        ----------
        item : int or slice or np.ndarray
            The index of the pair of atoms to get.

        Returns
        -------
        NeighborPair
            A new NeighborPair object containing the pair of atoms at the specified index.
        """
        if isinstance(item, int):
            return self.__class__(
                self._atoms, self._pairs[item].reshape(1, 2), self._distances[item]
            )

        return self.__class__(self._atoms, self.pairs[item], self.distances[item])

    def __contains__(self, other) -> bool:
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
        from ..utils.array_utils import issubset

        return issubset(other.pairs, self.pairs)

    def __add__(self, other):
        # TODO:
        # currently this is different from MDAnalysis.
        # The question to answer is if neighbor pairs should be unique or not.
        return self.union(other)

    def __sub__(self, other):
        return self.difference(other)

    def __or__(self, other):
        return self.symmetric_difference(other)

    def __eq__(self, other):
        return np.array_equal(self.pairs, other.pairs)

    def __and__(self, other):
        return self.intersection(other)

    def __xor__(self, other):
        return self.symmetric_difference(other)

    def __len__(self):
        return self.pairs.shape[0]

    def __str__(self) -> str:
        return f"<Lahuta NeighborPairs class containing {self.indices.size} atoms and {self.pairs.shape[0]} pairs>"

    def __repr__(self) -> str:
        return self.__str__()


class HBondNeighborPairs(NeighborPairs):
    """An extension of the `NeighborPairs` class for hydrogen bonds."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.hbond_array = find_hydrogen_bonded_atoms(self._atoms.universe.mol)
        self.hbond_angles = (
            None if kwargs.get("angles") is None else kwargs.get("angles")
        )
        self.result_array = (
            None if kwargs.get("result_array") is None else kwargs.get("result_array")
        )

    def hbond_distance_filter(
        self, col: int = 0, vdw_comp_factor: float = 0.1
    ) -> "HBondNeighborPairs":
        """Filter the pairs based on the distance between the hydrogen bonded atoms.

        Parameters
        ----------
        col : int
            The column of the hydrogen bonded atom indices in the `hbond_array`.
        vdw_comp_factor : float
            The van der Waals complementarity factor.

        Returns
        -------
        pairs : NeighborPairs
            The filtered neighbor pairs.
        """
        col2 = 2
        if col == 1:
            col2 = 1

        attr_col = getattr(self, f"col{col+1}")
        hbound_attr_col = getattr(self, f"col{col2}")

        tt_vdw_dist = attr_col.atoms.vdw_radii + VDW_RADII["H"] + vdw_comp_factor

        hbond_dist = self._get_hbond_distances(attr_col, hbound_attr_col)

        distance_condition = hbond_dist <= tt_vdw_dist[:, np.newaxis]
        tt_hbond_dist_pairs = self.pairs[np.any(distance_condition, axis=1)]
        tt_hbond_distances = self.distances[np.any(distance_condition, axis=1)]

        return self.__class__(
            self._atoms,
            tt_hbond_dist_pairs,
            tt_hbond_distances,
            angles=self.hbond_angles,
            result_array=self.result_array,
        )

    def hbond_angle_filter(
        self, col: int = 0, weak: bool = False
    ) -> "HBondNeighborPairs":
        """Filter the pairs based on the angle between the hydrogen bonded atoms.

        Parameters
        ----------
        col : int
            The column of the hydrogen bonded atom indices in the `hbond_array`.
        angle : float
            The angle in degrees.

        Returns
        -------
        pairs : NeighborPairs1
            The filtered neighbor pairs.
        """

        contact_type = "weak hbond" if weak else "hbond"
        col2 = 2
        if col == 1:
            col2 = 1

        attr_col = getattr(self, f"col{col+1}")
        hbound_attr_col = getattr(self, f"col{col2}")

        self.hbond_angles = self._get_hbond_angles(attr_col, hbound_attr_col)
        # self.hbond_angles = self._get_hbond_angles(hbound_attr_col, attr_col)

        idx = np.any(self.hbond_angles >= CONTACTS[contact_type]["angle rad"], axis=1)
        self._pairs = self._pairs[idx]
        self._distances = self._distances[idx]

        # self.result_array = np.full((len(self.hbond_angles), 4), np.nan)
        # self.result_array[:, :2] = self._pairs
        # self.result_array[:, 2] = np.nanmin(self.hbond_angles, axis=1)
        # self.result_array[:, 3] = self.distances

        return self.__class__(
            self._atoms,
            self._pairs,
            self.distances,
            angles=self.hbond_angles,
            result_array=self.result_array,
        )

    def _get_hbond_distances(self, attr_col, hbound_attr_col):
        """Get the distances between the hydrogen bonded atoms.

        Parameters
        ----------
        col : int
            The column index of the hydrogen bonded atom.

        Returns
        -------
        distances : np.ndarray
            The distances between the hydrogen bonded atoms.
        """

        col_atom_pos = attr_col.atoms.positions

        hbound_atom_indices = self.hbond_array[hbound_attr_col.atoms.indices]
        hbound_atom_pos = self._atoms.positions[hbound_atom_indices]

        hbound_atom_pos[hbound_atom_indices == 0] = np.nan
        distance_array = np.linalg.norm(
            col_atom_pos[:, np.newaxis, :] - hbound_atom_pos, axis=-1
        )

        return distance_array

    def _get_hbond_angles(self, col1, col2):
        atom1_pos = col1.atoms.positions
        atom2_pos = col2.atoms.positions

        hbound_atom_indices = self.hbond_array[col1.atoms.indices]
        hbound_atom_pos = self._atoms.positions[hbound_atom_indices]

        hbound_atom_pos[hbound_atom_indices == 0] = np.nan

        point_a = atom1_pos[:, np.newaxis, :]
        point_b = hbound_atom_pos
        point_c = atom2_pos[:, np.newaxis, :]

        return calculate_angle(point_a, point_b, point_c, degrees=False)
