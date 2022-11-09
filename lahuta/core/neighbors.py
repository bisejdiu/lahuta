"""
Placeholder for the neighbors module.
"""

from abc import abstractmethod
from typing import List, Union, Any
from typing_extensions import Protocol, Literal

import numpy as np

from ..utils.array import calculate_angle
from ..utils.atom_types import find_hydrogen_bonded_atoms
from ..config.defaults import CONTACTS, VDW_RADII

# from .groups import AtomGroup


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
        self.col1, self.col2 = uniatom[pairs[:, 0]], uniatom[pairs[:, 1]]

        # put kwargs in the object
        for key, value in kwargs.items():
            setattr(self, key, value)

        assert (
            self.col1.indices.size == self.col2.indices.size == self._distances.size
        ), "The `NeighborPairs` must have the number of indices."

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
        # return the parent class
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

    def hbond_distance_filter(self, col: int = 0, vdw_comp_factor: float = 0.1):
        raise NotImplementedError(
            "This method is not implemented for `NeighborPairs`. Please use `HBondNeighborPairs`."
        )

    @property
    def pairs(self) -> np.ndarray:
        """Get the pairs of atoms that are neighbors."""
        return self._pairs
        # return np.column_stack((self.col1.indices, self.col2.indices))

    @property
    def distances(self) -> np.ndarray:
        """Get the distances between the pairs of atoms that are neighbors."""
        return self._distances

    @property
    def indices(self) -> np.ndarray:
        """Get the indices of the atoms that are neighbors."""
        return np.unique([self.col1.indices, self.col2.indices])

    def __str__(self) -> str:
        return f"<Lahuta NeighborPairs class containing \
        {self.indices.size} atoms and {self.pairs.shape[0]} pairs>"

    def __repr__(self) -> str:
        return self.__str__()


class HBondNeighborPairs(NeighborPairs):
    """An extension of the `NeighborPairs` class for hydrogen bonds."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.hbond_array = find_hydrogen_bonded_atoms(self._atoms.universe.mol)
        self.angles = None if kwargs.get("angles") is None else kwargs.get("angles")
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
        # hbond_dist = self._get_hbond_distances(hbound_attr_col, attr_col)

        distance_condition = hbond_dist <= tt_vdw_dist[:, np.newaxis]
        tt_hbond_dist_pairs = self.pairs[np.any(distance_condition, axis=1)]
        tt_hbond_distances = self.distances[np.any(distance_condition, axis=1)]

        return self.__class__(
            self._atoms,
            tt_hbond_dist_pairs,
            tt_hbond_distances,
            angles=self.angles,
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

        self.angles = self._get_hbond_angles(attr_col, hbound_attr_col)
        # self.angles = self._get_hbond_angles(hbound_attr_col, attr_col)

        idx = np.any(self.angles >= CONTACTS[contact_type]["angle rad"], axis=1)
        self._pairs = self._pairs[idx]
        self._distances = self._distances[idx]

        # self.result_array = np.full((len(self.angles), 4), np.nan)
        # self.result_array[:, :2] = self._pairs
        # self.result_array[:, 2] = np.nanmin(self.angles, axis=1)
        # self.result_array[:, 3] = self.distances

        return self.__class__(
            self._atoms,
            self._pairs,
            self.distances,
            angles=self.angles,
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
