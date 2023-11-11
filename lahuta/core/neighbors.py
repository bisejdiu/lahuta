"""The primary class in this module is `NeighborPairs`, which represents pairs of atoms that are
considered as "neighbors" based on a certain distance criterion. This class provides a suite of
methods to manipulate, analyze, and export these pairs.
"""

from typing import TYPE_CHECKING, Any, Literal, Optional

import numpy as np
from numpy.typing import NDArray
from typing_extensions import Self

from lahuta.config.smarts import AVAILABLE_ATOM_TYPES
from lahuta.core.builder import AtomMapper, LabeledNeighborPairsBuilder
from lahuta.core.index_finder import IndexFinder
from lahuta.utils import array_utils as au

if TYPE_CHECKING:
    from Bio.Seq import Seq
    from pandas import DataFrame

    from lahuta._types.mdanalysis import AtomGroupType
    from lahuta.core.labeled_neighbors import LabeledNeighborPairs
    from lahuta.core.luni import Luni


__all__ = ["NeighborPairs"]


class NeighborPairs:
    """A class that manages pairs of atoms that are considered "neighbors" within a defined distance threshold.

    The `NeighborPairs` class stores atom pairs (identified by their indices) and the corresponding distances
    between them, which represent the concept of "neighbors" in the context of molecular simulations or structural
    biology. The class provides various functionalities to manipulate and interrogate these pairs, including
    methods to retrieve specific pairs, add or modify pairs, convert the data to different formats (e.g., pandas
    DataFrame or dictionary), and compute derived properties such as distances and indices.

    The class also implements several magic methods to support common operations like indexing, testing for
    membership, set-like operations (e.g., union, intersection), and equality testing. Furthermore, the class is
    designed to be extensible and supports the addition of custom annotations to the pairs.

    Args:
        luni (Luni): A Luni object containing the structure.

    Attributes:
        pairs (NDArray[np.int32]): A 2D numpy array of pairs of atom indices that are neighbors.
        distances (NDArray[np.float32]): A 1D numpy array of distances between the pairs of atoms.
        annotations (dict[str, NDArray[Any]]): A dictionary containing the annotations of the NeighborPairs object.
        partner1 (AtomGroupType): The first column of the pairs of atoms.
        partner2 (AtomGroupType): The second column of the pairs of atoms.
        indices (NDArray[np.int32]): A 2D numpy array of the indices of the pairs of atoms.

    Examples:
        ``` py
        >>> from lahuta import Luni, NeighborPairs
        >>> luni = Luni("path/to/structure.pdb")
        >>> np = NeighborPairs(luni)
        >>> np.set_neighbors(pairs, distances)
        >>> np.pairs
        [[1, 2], [3, 4]]
        >>> np.distances
        [1.0, 2.0]
        >>> np.annotations
        {}
        >>> np.partner1
        <AtomGroup with 2 atoms>
        >>> np.partner2
        <AtomGroup with 2 atoms>
        >>> np.indices
        [[1, 2], [3, 4]]
        ```
    """

    def __init__(self, luni: "Luni"):
        """Initialize the NeighborPairs object."""
        self.luni = luni
        self.atoms = luni.to("mda").universe.atoms

        self._pairs = np.array([], dtype=np.int32).reshape(0, 2)
        self._distances = np.array([], dtype=np.float32)
        self._annotations: dict[str, NDArray[Any]] = {}

    def _validate_inputs(self, pairs: NDArray[np.int32], distances: NDArray[np.float32]) -> None:
        """Validate that the provided pairs and distances arrays have the same first dimension.

        This internal method asserts that the first dimension of the `pairs` and `distances` arrays are equal.
        This is necessary to ensure that there is a one-to-one correspondence between the pairs of atoms
        and their respective distances. If the dimensions are not equal, the method raises an assertion error.

        Args:
            pairs (NDArray[np.int32]): An array containing the pairs of atoms.
            distances (NDArray[np.float32]): An array containing the distances between each pair of atoms.

        Raises:
            AssertionError: If the first dimension of the `pairs` and `distances` arrays are not equal.

        """
        message = (
            "The number of pairs and distances must be the same."
            f"Got {pairs.shape[0]} pairs and {distances.shape[0]} distances."
        )
        assert pairs.shape[0] == distances.shape[0], message

    @staticmethod
    def sort_inputs(
        pairs: NDArray[np.int32], distances: NDArray[np.float32]
    ) -> tuple[NDArray[np.int32], NDArray[np.float32]]:
        """Sorts the provided pairs and distances arrays based on the first column of the sorted pairs array.

        This method first sorts the `pairs` array along axis 1 and then sorts the `pairs` and `distances`
        arrays based on the first column of the sorted pairs array. This is done to ensure that the pairs and
        distances arrays are always in the same order, which is necessary for correctly associating
        each pair of atoms with its corresponding distance.

        Args:
            pairs (NDArray[np.int32]): An array containing the pairs of atoms.
            distances (NDArray[np.float32]): An array containing the distances between each pair of atoms.

        Returns:
            tuple[NDArray[np.int32], NDArray[np.float32]]: A tuple containing the sorted pairs and distances arrays.

        Example:
            ``` py
            pairs = np.array([[2, 1], [4, 3]])
            distances = np.array([1.0, 2.0])
            sorted_pairs, sorted_distances = NeighborPairs.sort_inputs(pairs, distances)
            print(sorted_pairs)
            print(sorted_distances)
            [[1, 2], [3, 4]]
            [2.0, 1.0]
            ```
        """
        pairs = np.sort(pairs, axis=1)
        indices = au.sorting_indices(pairs)

        return pairs[indices], distances[indices]

    def _get_pair_column(self, partner: int) -> "AtomGroupType":
        """Return the column of the pair of atoms depending on the value of partner."""
        return self.atoms[self.pairs[:, partner - 1]]

    def get_partners(self, partner: int) -> tuple["AtomGroupType", "AtomGroupType"]:
        """Return the columns of the pair of atoms depending on the value of partner."""
        partner2 = 1 if partner == 2 else 2

        return self._get_pair_column(partner), self._get_pair_column(partner2)

    def type_filter(self, atom_type: str, partner: int) -> Self:
        """Filter pairs based on atom types.

        The method selects pairs from the NeighborPairs object where the atoms have the specified type.
        The `partner` parameter specifies the column (1 or 2) from which the atom types are selected.

        Args:
            atom_type (str): Specifies the atom type. Must be one of:

                - 'carbonyl_oxygen'
                - 'weak_hbond_donor'
                - 'pos_ionisable'
                - 'carbonyl_carbon'
                - 'hbond_acceptor'
                - 'hbond_donor'
                - 'neg_ionisable'
                - 'weak_hbond_acceptor'
                - 'xbond_acceptor'
                - 'aromatic'
                - 'hydrophobe'

            partner (int): The column for atom type selection. Can be either 1 or 2.

        Returns:
            A NeighborPairs object containing the pairs that meet the atom type filter.
        """
        atom_type_col_num = AVAILABLE_ATOM_TYPES[atom_type.upper()]
        nonzeros: NDArray[np.int32] = self.luni.atom_types.getcol(atom_type_col_num).nonzero()[0]
        mask = np.in1d(self.pairs[:, partner - 1], nonzeros)

        return self.new(self.pairs[mask], self.distances[mask])

    def index_filter(
        self,
        indices: NDArray[np.int32],
        partner: int,
    ) -> Self:
        """Select pairs based on the atom indices.

        The method selects pairs from the NeighborPairs object where the atoms have the specified indices.
        The `partner` parameter specifies the column (1 or 2) from which the atom indices are selected.

        Args:
            indices (NDArray[np.int32]): The atom indices to select.
            partner (int): The column to select the atom indices from. It can be either 1 or 2.

        Returns:
            A NeighborPairs object containing the pairs that meet the index filter.
        """
        mask = np.in1d(self.pairs[:, partner - 1], indices)
        return self.new(self.pairs[mask], self.distances[mask])

    def distance_filter(self, distance: float) -> Self:
        """Select pairs based on the distance.

        The method selects pairs from the NeighborPairs object where the distances between the atoms are
        less than or equal to the specified distance.

        Args:
            distance (float): The distance to select.

        Returns:
            A NeighborPairs object containing the pairs that meet the distance filter.
        """
        mask = self.distances <= distance
        return self.new(self.pairs[mask], self.distances[mask])

    def numeric_filter(self, array: NDArray[np.float32], cutoff: float) -> Self:
        """Select pairs based on a numeric cutoff.

        The method selects pairs from the NeighborPairs object where the values in the specified array are less than or
        equal to the cutoff (if `lte` is True) or greater than the cutoff (if `lte` is False).

        Args:
            array (NDArray[np.float32]): The array containing the values to compare with the cutoff.
            cutoff (float): The cutoff value for the filter.

        Returns:
            A NeighborPairs object containing the pairs that meet the numeric filter.
        """
        mask = array <= cutoff
        return self.new(self.pairs[mask], self.distances[mask])

    def radius_filter(self, radius: float, partner: int) -> Self:
        """Select pairs based on the radius.

        The method selects pairs from the NeighborPairs object where the van der Waals radii of the atoms
        are less than or equal to the specified radius. The `partner` parameter specifies the column (1 or 2)
        from which the radii are selected.

        Args:
            radius (float): The radius to select.
            partner (int): The column to select the radii from. It can be either 1 or 2.

        Returns:
            A NeighborPairs object containing the pairs that meet the radius filter.
        """
        col_func = self._get_pair_column(partner)
        mask = col_func.atoms.vdw_radii <= radius

        return self.new(self.pairs[mask], self.distances[mask])

    def map(self, seq: "Seq") -> "LabeledNeighborPairs":
        """Map the `pairs` indices to indices in the multiple sequence alignment.

        The method maps the indices in the `pairs` array to indices in the multiple sequence alignment
        using the specified sequence ID.

        Args:
            seq (Bio.Seq): The sequence to map the indices to.

        Returns:
            A NeighborPairs object containing the mapped pairs.
        """
        atom_mapper = AtomMapper(self.atoms)
        builder = LabeledNeighborPairsBuilder(atom_mapper)
        return builder.build(self.pairs, seq)

    def backmap(self, seq: "Seq", pairs: NDArray[np.void]) -> Self:
        """Map the `pairs` indices to indices in the structure.

        The method maps the indices in the `pairs` array to indices in the structure
        using the specified sequence ID.

        Args:
            seq (Bio.Seq): The sequence to map the indices to.
            pairs (NDArray[np.void]): The mapped pairs to backmap.

        Returns:
            A NeighborPairs object containing the backmapped pairs.
        """
        atom_mapper = AtomMapper(self.atoms)
        mapped_pairs = LabeledNeighborPairsBuilder(atom_mapper).build(self.pairs, seq).pairs
        index_finder = IndexFinder(mapped_pairs)
        mask = index_finder.find_indices(pairs)
        return self.new(self.pairs[mask], self.distances[mask])

    def intersection(self, other: Self) -> Self:
        """Return the intersection of two NeighborPairs objects.

        The method calculates the intersection of the pairs from `self` and `other`, and then returns a new
        NeighborPairs object that contains the intersecting pairs along with their corresponding distances.

        Args:
            other: The other NeighborPairs object.

        Returns:
            intersected_pairs: A NeighborPairs object containing the pairs and their corresponding distances
                            that are common between `self` and `other`.

        Example:
            ``` py
            >>> np1 = NeighborPairs(...)
            >>> np2 = NeighborPairs(...)
            # 2 ways to get the intersection
            >>> np_intersected = np1 & np2
            >>> np_intersected = np1.intersection(np2)
            ```
        """
        mask = au.intersection(self.pairs, other.pairs)
        return self.new(self.pairs[mask], self.distances[mask])

    def union(self, other: Self) -> Self:
        """Return the union of two NeighborPairs objects.

        The method finds the union of the pairs from `self` and `other`. It also ensures that the distances
        in the resulting object correspond to the union pairs.

        Args:
            other: The other NeighborPairs object to be unified with.

        Returns:
            pairs: A NeighborPairs object containing the union of the pairs from `self` and `other`, and
                with corresponding distances.

        Example:
            ``` py
            >>> np1 = NeighborPairs(...)
            >>> np2 = NeighborPairs(...)
            # 2 ways to get the union
            >>> np_union = np1 + np2
            >>> np_union = np1.union(np2)
            ```
        """
        pairs, indices = au.union(self.pairs, other.pairs)
        distances = np.concatenate((self.distances, other.distances), axis=0)[indices]

        return self.new(pairs, distances)

    def difference(self, other: Self) -> Self:
        """Return the difference between two NeighborPairs objects.

        The method calculates the difference between the pairs from `self` and `other`, then returns a new
        NeighborPairs object that contains the pairs from `self` that are not in `other`, along with their
        corresponding distances.

        Args:
            other: The other NeighborPairs object.

        Returns:
            difference_pairs: A NeighborPairs object containing the pairs and their corresponding distances
                            from `self` that are not in `other`.

        Example:
            ``` py
            >>> np1 = NeighborPairs(...)
            >>> np2 = NeighborPairs(...)
            # 2 ways to get the difference
            >>> np_diff = np1 - np2
            >>> np_diff = np1.difference(np2)
            ```

        """
        mask = au.difference(self.pairs, other.pairs)

        return self.new(self.pairs[mask], self.distances[mask])

    def symmetric_difference(self, other: Self) -> Self:
        """Return the symmetric difference of two NeighborPairs objects.

        This method creates a new `NeighborPairs` object that contains pairs and distances
        that are unique to `self` or `other`, but not both.

        Args:
            other: The other NeighborPairs object.

        Returns:
            A NeighborPairs object containing the symmetric difference of the two NeighborPairs objects.

        Example:
            ``` py
            >>> np1 = NeighborPairs(...)
            >>> np2 = NeighborPairs(...)
            # 2 ways to get the symmetric difference
            >>> np_sym_diff = np1 | np2
            >>> np_sym_diff = np1.symmetric_difference(np2)
            ```
        """
        mask_a, mask_b = au.symmetric_difference(self.pairs, other.pairs)

        pairs = np.concatenate((self.pairs[mask_a], other.pairs[mask_b]), axis=0)
        distances = np.concatenate((self.distances[mask_a], other.distances[mask_b]), axis=0)

        return self.new(pairs, distances)

    def isdisjoint(self, other: Self) -> bool:
        """Check if the intersection of two NeighborPairs objects is null.

        This method checks whether the intersection of the two NeighborPairs objects is null,
        thus determining if the two objects are disjoint.

        Args:
            other (NeighborPairs): The other NeighborPairs object.

        Returns:
            bool: True if the two NeighborPairs objects are disjoint
                    (i.e., have no common pairs), and False otherwise.

        Example:
            ``` py
            >>> np1 = NeighborPairs(...)
            >>> np2 = NeighborPairs(...)
            >>> np1.isdisjoint(np2)
            True
            ```
        """
        return au.isdisjoint(self.pairs, other.pairs)

    def issubset(self, other: Self) -> bool:
        """Check if all elements (pairs) of a NeighborPairs object are found in another NeighborPairs object.

        This method checks whether every pair of atoms from the current NeighborPairs object
        is also present in the other NeighborPairs object, thus determining if this object is a subset of the 'other'.

        Args:
            other (NeighborPairs): The other NeighborPairs object.

        Returns:
            bool: True if every pair in the current NeighborPairs object is found
                    in the other NeighborPairs object, and False otherwise.

        Example:
            ``` py
            >>> np1 = NeighborPairs(...)
            >>> np2 = NeighborPairs(...)
            >>> np1.issubset(np2)
            ```
        """
        return au.issubset(self.pairs, other.pairs)

    def issuperset(self, other: Self) -> bool:
        """Determine if all pairs from another NeighborPairs object are found in this object.

        This method checks whether every pair of atoms from the 'other' NeighborPairs object
        is also present in this object, thus determining if this object is a superset of the 'other'.

        Args:
            other (NeighborPairs): The other NeighborPairs object.

        Returns:
            bool: True if all pairs from 'other' are found in this object, False otherwise.

        Example:
            ``` py
            >>> np1 = NeighborPairs(...)
            >>> np2 = NeighborPairs(...)
            >>> np1.issuperset(np2)
            True
            ```
        """
        return au.issuperset(self.pairs, other.pairs)

    def isequal(self, other: Self) -> bool:
        """Check if this NeighborPairs object is equal to another.

        Two NeighborPairs objects are considered equal if they contain exactly the same pairs.

        Args:
            other (NeighborPairs): The other NeighborPairs object.

        Returns:
            bool: True if the two NeighborPairs objects contain the same pairs, False otherwise.

        Example:
            ``` py
            >>> np1 = NeighborPairs(...)
            >>> np2 = NeighborPairs(...)
            >>> np1.isequal(np2)
            True
            ```
        """
        return au.isequal(self.pairs, other.pairs)

    def isunique(self) -> bool:
        """Check if all pairs in this NeighborPairs object are unique.

        This method checks if all pairs in this NeighborPairs object are unique, i.e., there are no duplicate pairs.

        Returns:
            bool: True if all pairs in this object are unique, False otherwise.

        Example:
            ``` py
            >>> np = NeighborPairs(...)
            >>> np.isunique()
            False
            ```
        """
        return au.isunique(self.pairs)

    def is_strict_subset(self, other: Self) -> bool:
        """Check if all pairs of this NeighborPairs object are in another, and the two sets are not equal.

        A strict subset has all pairs in the 'other' object but the two sets are not identical.

        Args:
            other (NeighborPairs): The other NeighborPairs object.

        Returns:
            bool: True if this object is a strict subset of 'other', False otherwise.

        Example:
            ``` py
            >>> np1 = NeighborPairs(...)
            >>> np2 = NeighborPairs(...)
            >>> np1.is_strict_subset(np2)
            True
            ```
        """
        return au.is_strict_subset(self.pairs, other.pairs)

    def is_strict_superset(self, other: Self) -> bool:
        """Check if all pairs of another NeighborPairs object are in this one, and the two sets are not equal.

        A strict superset has all pairs from the 'other' object but the two sets are not identical.

        Args:
            other (NeighborPairs): The other NeighborPairs object.

        Returns:
            bool: True if this object is a strict superset of 'other', False otherwise.

        Example:
            ``` py
            >>> np1 = NeighborPairs(...)
            >>> np2 = NeighborPairs(...)
            >>> np1.is_strict_superset(np2)
            True
            ```
        """
        return au.is_strict_superset(self.pairs, other.pairs)

    def new(self, pairs: NDArray[np.int32], distances: NDArray[np.float32]) -> Self:
        """Return a new NeighborPairs object that is a copy of the current object,
        but with specified pairs and distances.

        Args:
            pairs (NDArray[np.int32]): The atom pairs for the new object.
            distances (NDArray[np.float32]): The corresponding distances for the new object.

        Returns:
            A new NeighborPairs object with the provided pairs and distances.
        """
        cls = type(self)
        new = cls.__new__(cls)
        new.luni = self.luni
        new.atoms = self.atoms
        new.set_neighbors(pairs, distances, sort=False)
        new.annotations = {}  # reset annotations
        return new

    def plot(self, which: Literal["matching", "full"] = "matching", half_only: bool = False) -> None:
        """Plot the contact map of the NeighborPairs object.

        Args:
            which (str, optional): Which contact map to plot ('matching' or 'full'). Defaults to 'matching'.
            half_only (bool, optional): Whether to plot only the upper half of the contact map. Defaults to False.
        """
        from lahuta.viz.contact_matrix import ContactMap
        return ContactMap(self.pairs).plot(which, half_only)

    @property
    def annotations(self) -> dict[str, NDArray[Any]]:
        """Get the annotations of the NeighborPairs object.

        Returns:
            A dictionary containing the annotations of the NeighborPairs object.
        """
        return self._annotations

    @annotations.setter
    def annotations(self, annotations: dict[str, NDArray[Any]]) -> None:
        """Set the annotations of the NeighborPairs object.

        Args:
            annotations (dict[str, NDArray[Any]]): A dictionary containing the annotations to be set.
        """
        self._annotations = annotations

    def add_annotations(self, annotations: dict[str, NDArray[Any]]) -> None:
        """Add annotations to the existing NeighborPairs object.

        Args:
            annotations (dict[str, NDArray[Any]]): A dictionary containing the annotations to be added.
        """
        for value in annotations.values():
            assert len(value) == self.pairs.shape[0]

        self._annotations.update(annotations)

    def to_frame(
        self,
        df_format: Literal["compact", "expanded"] = "expanded",
        annotations: bool = False,
    ) -> "DataFrame":
        """Convert the NeighborPairs object to a pandas DataFrame.

        The method provides two formatting options. The 'compact' format contains two columns
        for atom indices and one column for distances. The 'expanded' format contains four columns
        for atom indices (two columns for each atom pair) and one column for distances.
        If `annotations` is True, the resulting DataFrame will also include annotation columns.

        Args:
            df_format (str, optional): The format of the DataFrame. It can be either "compact" or "expanded".
                                        Defaults to "expanded".
            annotations (bool, optional): Whether to include annotations in the DataFrame. Defaults to False.

        Returns:
            A pandas DataFrame containing the atom pairs and their distances.
        """
        if annotations:
            return self._create_df(df_format, self.annotations)

        return self._create_df(df_format)

    def to_dict(self, df_format: Literal["compact", "expanded"] = "expanded") -> dict[str, Any]:
        """Convert the NeighborPairs object to a dictionary.

        The method converts the NeighborPairs object to a pandas DataFrame
        (using the specified format), and then converts that DataFrame to a dictionary.

        Args:
            df_format (str, optional): The format of the DataFrame. It can be either "compact" or "expanded".
                                        Defaults to "expanded".

        Returns:
            A dictionary representation of the NeighborPairs object.
        """
        return self._create_df(df_format).to_dict(orient="list")  # type: ignore

    def _create_df(
        self,
        df_format: Literal["compact", "expanded"] = "expanded",
        annotations: Optional[dict[str, NDArray[Any]]] = None,
    ) -> "DataFrame":
        """Create a pandas DataFrame from the NeighborPairs object.

        Args:
            df_format (str, optional): The format of the DataFrame.
            It can be either "compact" or "expanded". Defaults to "expanded".
            annotations: Dictionary containing additional information
                        to be added to DataFrame. Defaults to None.

        Returns:
            A pandas DataFrame representing the NeighborPairs object.
        """
        from lahuta.writers.frame_writer import DataFrameWriter
        return DataFrameWriter(self, df_format, annotations).create()

    def vmd_exporter(self, sphere_resolution: int = 20, save_to_file: bool = False) -> Optional[str]:
        """Export a TCL script to visualize the neighbor pairs in VMD.

        Args:
            sphere_resolution (int, optional): The resolution of the spheres. Defaults to 20.
            save_to_file (bool, optional): Flag to save the script to a file. Defaults to False.

        Returns:
            str | None: The TCL script to visualize the neighbor pairs in VMD.
        """
        from lahuta.writers.exporters import VMDExporter
        exporter = VMDExporter(self.pairs)
        return exporter.export(sphere_resolution=sphere_resolution, save_to_file=save_to_file)

    def _neighborpairs_equal(self, other: Self) -> bool:
        """Check if another NeighborPairs object is equal to this one.

        Two NeighborPairs objects are considered equal if they have the same atom pairs
        and the same corresponding distances.

        Args:
            other (NeighborPairs): The other NeighborPairs object to compare to this one.

        Returns:
            True if the other object is equal to this one; False otherwise.
        """
        indices1: NDArray[np.int32] = np.lexsort((self.pairs[:, 1], self.pairs[:, 0]))
        indices2: NDArray[np.int32] = np.lexsort((other.pairs[:, 1], other.pairs[:, 0]))

        # Sort each array using the indices
        pairs, dists = self.pairs[indices1], self.distances[indices1]
        other_pairs, other_dists = other.pairs[indices2], other.distances[indices2]

        return np.array_equal(pairs, other_pairs) and np.allclose(dists, other_dists)

    @property
    def partner1(self) -> "AtomGroupType":
        """Get the first partner of the pairs of indices of atoms that are neighbors.

        Returns:
            The first partner of the atom pairs.
        """
        return self._get_pair_column(1)

    @property
    def partner2(self) -> "AtomGroupType":
        """Get the second partner of the pairs of indices of atoms that are neighbors.

        Returns:
            The second partner of the atom pairs.
        """
        return self._get_pair_column(2)

    @property
    def pairs(self) -> NDArray[np.int32]:
        """Get the pairs of atoms that are neighbors.

        Returns:
            An array containing the pairs of indices of neighboring atoms.
        """
        return self._pairs

    @property
    def neighbors(self) -> tuple[NDArray[np.int32], NDArray[np.float32]]:
        """Get the neighbors.

        Returns:
            An array containing the indices of the neighbors.
        """
        return self.pairs, self.distances

    @neighbors.setter
    def neighbors(self, pairs: NDArray[np.int32], distances: NDArray[np.float32]) -> None:
        """Set the neighbors.

        Args:
            pairs (NDArray[np.int32]): An array containing the pairs of atoms.
            distances (NDArray[np.float32]): An array containing the distances between each pair of atoms.
        """
        self.set_neighbors(pairs, distances)

    def get_neighbors(self) -> tuple[NDArray[np.int32], NDArray[np.float32]]:
        """Get the neighbors.

        Returns:
            An array containing the indices of the neighbors.
        """
        return self.pairs, self.distances

    def set_neighbors(self, pairs: NDArray[np.int32], distances: NDArray[np.float32], sort: bool = True) -> None:
        """Set the neighbors.

        Args:
            pairs (NDArray[np.int32]): An array containing the pairs of atoms.
            distances (NDArray[np.float32]): An array containing the distances between each pair of atoms.
            sort (bool, optional): Whether to sort the pairs and distances arrays. Defaults to True.
        """
        self._validate_inputs(pairs, distances)
        if sort:
            self._pairs, self._distances = NeighborPairs.sort_inputs(pairs, distances)
        else:
            self._pairs, self._distances = pairs, distances

    @property
    def distances(self) -> NDArray[np.float32]:
        """Get the distances between the pairs of indices of atoms that are neighbors.

        Returns
            An array containing the distances between the pairs of indices of neighboring atoms.
        """
        return self._distances

    @property
    def names(self) -> NDArray[np.str_]:
        """Get the names of the atoms that are neighbors.

        Returns
            An array containing the atom names of the neighboring atoms.
        """
        return self.atoms[self.pairs].names

    @property
    def resnames(self) -> NDArray[np.str_]:
        """Get the residue names of the atoms that are neighbors.

        Returns
            An array containing the residue names of the neighboring atoms.
        """
        return self.atoms[self.pairs].resnames

    @property
    def resids(self) -> NDArray[np.int32]:
        """Get the residue IDs of the atoms that are neighbors.

        Returns
            An array containing the residue IDs of the neighboring atoms.
        """
        return self.atoms[self.pairs].resids

    @property
    def chainIDs(self) -> NDArray[np.int32]:
        """Get the chain IDs of the atoms that are neighbors.

        Returns
            An array containing the chain IDs of the neighboring atoms.
        """
        return self.atoms[self.pairs].chainIDs

    @property
    def indices(self) -> NDArray[np.int32]:
        """Get the indices of the atoms that are neighbors.

        Returns
            An array containing the indices of the neighboring atoms.
        """
        return self.atoms[self.pairs].indices

    @property
    def labels(self) -> NDArray[np.void]:
        """Get the pairs of atoms that are neighbors.

        Returns:
            An array containing the pairs of indices of neighboring atoms.
        """
        struct_array = LabeledNeighborPairsBuilder.create_empty_struct_array(self.atoms.n_atoms)
        struct_array["names"] = self.atoms.names
        struct_array["resnames"] = self.atoms.resnames
        struct_array["resids"] = self.atoms.resids

        return struct_array[self.pairs]

    def __getitem__(self, item: int | slice | NDArray[np.int32]) -> Self:
        """Retrieve the neighbor pairs at the specified index or indices.

        This method allows accessing the neighbor pairs similar to elements in a list.
        For an integer input, it returns a NeighborPairs object containing a single pair,
        while for a slice or an array of indices, it returns a NeighborPairs object with the corresponding pairs.

        Args:
            item (int, slice, or ndarray): An integer, slice, or array of integers indicating
                                        the index/indices of the pair(s).

        Returns:
            NeighborPairs: A new NeighborPairs object containing the specified pair(s) and
                            their corresponding distance(s).
        """
        if isinstance(item, int):
            return self.new(
                self.pairs[item].reshape(-1, 2),
                np.array([self.distances[item]]),
            )

        return self.new(self.pairs[item], self.distances[item])

    def __contains__(self, other: Self) -> bool:
        """Check whether all pairs in the given NeighborPairs object are also present in this NeighborPairs object.

        This method allows using the Python built-in `in` keyword to check for the presence of pairs.
        This method add support for the `in` operator.

        Args:
            other (NeighborPairs): The NeighborPairs object to check.

        Returns:
            bool: True if all pairs from the 'other' NeighborPairs object are found in this one; False otherwise.

        Raises:
            NotImplemented: If the 'other' object is not an instance of the NeighborPairs class.
        """
        if other.__class__ != self.__class__:
            return NotImplemented

        return au.issubset(other.pairs, self.pairs)

    def __add__(self, other: Self) -> Self:
        """Combine this NeighborPairs object with another one.

        The resulting NeighborPairs object is the union of the two sets, containing all unique pairs from both.
        This method adds support for the `+` operator.

        Args:
            other (NeighborPairs): Another NeighborPairs object.

        Returns:
            NeighborPairs: A new NeighborPairs object that is the union of this one and the 'other'.

        Raises:
            NotImplemented: If the 'other' object is not an instance of the NeighborPairs class.
        """
        if other.__class__ != self.__class__:
            return NotImplemented
        return self.union(other)

    def __sub__(self, other: Self) -> Self:
        """Get the pairs in this NeighborPairs object that are not in the 'other'.

        The resulting NeighborPairs object is the difference of the two sets,
        containing pairs present in this object but not in the 'other'.
        This method adds support for the `-` operator.

        Args:
            other (NeighborPairs): Another NeighborPairs object.

        Returns:
            NeighborPairs: A new NeighborPairs object that is the difference of this one and the 'other'.

        Raises:
            NotImplemented: If the 'other' object is not an instance of the NeighborPairs class.
        """
        if other.__class__ != self.__class__:
            return NotImplemented
        return self.difference(other)

    def __or__(self, other: Self) -> Self:
        """Get the pairs that are in either this NeighborPairs object or the 'other', but not in both.

        The resulting NeighborPairs object is the symmetric difference of the two sets,
        containing pairs present in either this object or the 'other', but not in both.
        This method adds support for the `|` operator.

        Args:
            other (NeighborPairs): Another NeighborPairs object.

        Returns:
            NeighborPairs: A new NeighborPairs object that is the symmetric difference of this one and the 'other'.

        Raises:
            NotImplemented: If the 'other' object is not an instance of the NeighborPairs class.
        """
        if other.__class__ != self.__class__:
            return NotImplemented
        return self.symmetric_difference(other)

    def __eq__(self, other: Any) -> bool:  # noqa: PYI032, ANN401
        """Check if this NeighborPairs object is equal to the 'other'.

        Equality is based on the pairs and their distances.
        This method adds support for the `==` operator.

        Args:
            other (Any): Another object.

        Returns:
            bool: True if 'other' is an identical NeighborPairs object; False otherwise.

        Raises:
            NotImplemented: If the 'other' object is not an instance of the NeighborPairs class.
        """
        if other.__class__ != self.__class__:
            return NotImplemented
        return self._neighborpairs_equal(other)

    def __and__(self, other: Self) -> Self:
        """Get the pairs that are common to both this NeighborPairs object and the 'other'.

        The resulting NeighborPairs object is the intersection of the two sets,
        containing pairs present in both this object and the 'other'.
        This method adds support for the `&` operator.

        Args:
            other (NeighborPairs): Another NeighborPairs object.

        Returns:
            NeighborPairs: A new NeighborPairs object that is the intersection of this one and the 'other'.

        Raises:
            NotImplemented: If the 'other' object is not an instance of the NeighborPairs class.
        """
        if other.__class__ != self.__class__:
            return NotImplemented
        return self.intersection(other)

    def __xor__(self, other: Self) -> Self:
        """Get the pairs that are in either this NeighborPairs object or the 'other', but not in both.

        This method behaves similarly to the `__or__` method.
        This method adds support for the `^` operator.

        Args:
            other (NeighborPairs): Another NeighborPairs object.

        Returns:
            NeighborPairs: A new NeighborPairs object that is the symmetric difference of this one and the 'other'.

        Raises:
            NotImplemented: If the 'other' object is not an instance of the NeighborPairs class.
        """
        if other.__class__ != self.__class__:
            return NotImplemented
        return self.symmetric_difference(other)

    def __lt__(self, other: Self) -> bool:
        if other.__class__ != self.__class__:
            return NotImplemented
        return self.is_strict_subset(other)

    def __le__(self, other: Self) -> bool:
        if other.__class__ != self.__class__:
            return NotImplemented

        return self.issubset(other)

    def __gt__(self, other: Self) -> bool:
        if other.__class__ != self.__class__:
            return NotImplemented

        return self.is_strict_superset(other)

    def __ge__(self, other: Self) -> bool:
        if other.__class__ != self.__class__:
            return NotImplemented

        return self.issuperset(other)

    def __ne__(self, other: Any) -> bool:  # noqa: PYI032, ANN401
        if other.__class__ != self.__class__:
            return NotImplemented
        return not self._neighborpairs_equal(other)

    def __len__(self) -> int:
        """Get the number of pairs in this NeighborPairs object."""
        return self.pairs.shape[0]

    def __str__(self) -> str:
        unique_indices = np.unique(self.indices.ravel())
        return f"<Lahuta NeighborPairs class containing {unique_indices.size} atoms and {self.pairs.shape[0]} pairs>"

    def __repr__(self) -> str:
        return self.__str__()
