"""Defines and manages the concept of "neighbors" within a molecular or biological context.

The primary class in this module is `LabeledNeighborPairs`, which represents pairs of atoms that are 
considered as "neighbors" based on a certain distance criterion. This class provides a suite of 
methods to manipulate, analyze, and export these pairs.

"""
import warnings
from typing import Any, Callable, Optional, Type, TypeVar

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from typing_extensions import Self

from lahuta._types.mdanalysis import AtomGroupType
from lahuta.msa.encoder import encode_labels
from lahuta.utils import array_utils as au

__all__ = ["LabeledNeighborPairs"]

T = TypeVar("T")


class LabeledNeighborPairs:
    """Manages pairs of atoms that are considered "neighbors" within a defined distance threshold.

    The `LabeledNeighborPairs` class stores atom pairs (identified by their indices) and the corresponding distances
    between them, which represent the concept of "neighbors" in the context of molecular simulations or structural
    biology. The class provides various functionalities to manipulate and interrogate these pairs, including
    methods to retrieve specific pairs, add or modify pairs, convert the data to different formats (e.g., pandas
    DataFrame or dictionary), and compute derived properties such as distances and indices.

    The class also implements several magic methods to support common operations like indexing, testing for
    membership, set-like operations (e.g., union, intersection), and equality testing. Furthermore, the class is
    designed to be extensible and supports the addition of custom annotations to the pairs.

    Args:
        pairs (NDArray[np.void]): The array of pairs.
        mask1 (Optional[NDArray[np.bool_]], optional): A mask to select the first atom of each pair. Defaults to None.
        mask2 (Optional[NDArray[np.bool_]], optional): A mask to select the second atom of each pair. Defaults to None.


    Attributes:
        _pairs (NDArray[np.void]): The array of pairs.
        _mask1 (NDArray[np.bool_]): A mask to select the first atom of each pair.
        _mask2 (NDArray[np.bool_]): A mask to select the second atom of each pair.

    """

    def __init__(
        self,
        pairs: NDArray[np.void],
        mask1: Optional[NDArray[np.bool_]] = None,
        mask2: Optional[NDArray[np.bool_]] = None,
    ):
        upairs = np.unique(pairs, axis=0)
        self._pairs = pairs if pairs.shape[0] == upairs.shape[0] else upairs
        self._mask1 = mask1 if mask1 is not None else np.ones(self._pairs.shape[0], dtype=bool)
        self._mask2 = mask2 if mask2 is not None else np.ones(self._pairs.shape[0], dtype=bool)

    def type_filter(self, *_: T) -> Self:
        """Not implemented."""
        raise NotImplementedError(f"{self.__class__.__name__} does not support type filtering.")

    def index_filter(self, *_: T) -> Self:
        """Not implemented."""
        raise NotImplementedError(f"{self.__class__.__name__} does not support index filtering.")

    def distance_filter(self, *_: T) -> Self:
        """Not implemented."""
        raise NotImplementedError(f"{self.__class__.__name__} does not support distance filtering.")

    def numeric_filter(self, *_: T) -> Self:
        """Not implemented."""
        raise NotImplementedError(f"{self.__class__.__name__} does not support numeric filtering.")

    def radius_filter(self, *_: T) -> Self:
        """Not implemented."""
        raise NotImplementedError(f"{self.__class__.__name__} does not support radius filtering.")

    def hbond_distance_filter(self, *_: T) -> Self:
        """Not implemented."""
        raise NotImplementedError(f"{self.__class__.__name__} does not support hydrogen bond distance filtering.")

    def hbond_angle_filter(self, *_: T) -> Self:
        """Not implemented."""
        raise NotImplementedError(f"{self.__class__.__name__} does not support hydrogen bond angle filtering.")

    def _filter(self, mode: str, inverse: bool = False, **kwargs: list[str]) -> Self:
        masks: list[NDArray[np.bool_]] = []
        for mask in [self._mask1, self._mask2]:
            new_mask = np.copy(mask)
            assert self._pairs.dtype.names is not None
            for key, values in kwargs.items():
                if key not in self._pairs.dtype.names:
                    raise ValueError(f"Field {key} not found in pairs")
                field_mask = np.isin(self._pairs[key][:, 0 if mask is self._mask1 else 1], values)
                new_mask &= ~field_mask if mode == "exclude" or inverse else field_mask

            masks.append(new_mask)
        return self.__class__(self._pairs, *masks)

    def select(self, **kwargs: list[str]) -> Self:
        """Select pairs based on their attributes.

        Possible attributes are: `names`, `resnames`, and `resids`. The values of these attributes
        should be provided as lists of strings. The method returns a new LabeledNeighborPairs object
        containing the pairs that match the provided attributes.

        Args:
            **kwargs: The attributes to filter by. Possible attributes are: `names`, `resnames`, and `resids`.

        Returns:
            LabeledNeighborPairs: A new LabeledNeighborPairs object containing the pairs that \
                match the provided attributes.

        Example:
            ``` py
            >>> np = LabeledNeighborPairs(...)
            >>> np.select(names=['CA', 'CB']) # select all pairs with CA or CB for names
            >>> np.select(resnames=['ALA', 'GLY']) # select all pairs with ALA or GLY for resnames
            >>> np.select(names=['CA', 'CB'], resnames=['ALA', 'GLY'])
            ```
        """
        return self._filter("select", False, **kwargs)

    def exclude(self, **kwargs: list[str]) -> Self:
        """Exclude pairs based on their attributes.

        Possible attributes are: `names`, `resnames`, and `resids`. The values of these attributes
        should be provided as lists of strings. The method returns a new LabeledNeighborPairs object
        containing the pairs that do not match the provided attributes.

        Args:
            **kwargs: The attributes to filter by. Possible attributes are: `names`, `resnames`, and `resids`.

        Returns:
            LabeledNeighborPairs: A new LabeledNeighborPairs object containing the pairs that \
                do not match the provided attributes.

        Example:
            ``` py
            >>> np = LabeledNeighborPairs(...)
            >>> np.exclude(names=['CA', 'CB']) # exclude all pairs with CA or CB for names
            >>> np.exclude(resnames=['ALA', 'GLY']) # exclude all pairs with ALA or GLY for resnames
            >>> np.exclude(names=['CA', 'CB'], resnames=['ALA', 'GLY'])
            ```

        """
        return self._filter("exclude", False, **kwargs)

    def inverse(self) -> Self:
        """Invert the selection.

        This method returns a new LabeledNeighborPairs object containing the pairs that are not selected
        by the current object.

        Returns:
            LabeledNeighborPairs: A new LabeledNeighborPairs object containing the pairs \
                that are not selected by the current object.

        Example:
            ``` py
            >>> np = LabeledNeighborPairs(...)
            >>> np.select(names=['CA', 'CB']).inverse() # select all pairs that do not have CA or CB for names
            ```
        """
        return self.__class__(self._pairs, ~self._mask1, ~self._mask2)

    @staticmethod
    def remove_duplicates(pairs: NDArray[np.void]) -> NDArray[np.void]:
        """Remove duplicate pairs.

        This method removes duplicate pairs from the provided array of pairs. The method returns a new array
        containing the unique pairs.

        Args:
            pairs (NDArray[np.str_]): The array of pairs to remove duplicates from.

        Returns:
            NDArray[np.void]: A new array containing the unique pairs.

        """
        return np.unique(pairs, axis=0)

    def _apply_func(self, field: str, func: Callable[[NDArray[np.str_]], str | np.str_]) -> Self:
        assert self._pairs.dtype.names is not None
        if field not in self._pairs.dtype.names:
            raise ValueError(f"Field {field} not found in pairs")

        new_data = np.copy(self._pairs)
        new_data[field][self._mask1, 0] = func(new_data[field][self._mask1, 0])
        new_data[field][self._mask2, 1] = func(new_data[field][self._mask2, 1])

        return self.__class__(self.remove_duplicates(new_data))

    def remove(self, field: str) -> Self:
        """Remove an attribute from the pairs.

        This method removes the specified attribute from the pairs. The method returns a new LabeledNeighborPairs object
        with the specified attribute removed.

        Args:
            field (str): The attribute to remove.

        Returns:
            LabeledNeighborPairs: A new LabeledNeighborPairs object with the specified attribute removed.

        Example:
            ``` py
            >>> np = LabeledNeighborPairs(...)
            >>> np.remove('names')
            ```
        """
        return self._apply_func(field, lambda _: "")

    def rename(self, field: str, func: Callable[[NDArray[np.str_]], str | np.str_]) -> Self:
        """Rename an attribute of the pairs.

        This method renames the specified attribute of the pairs using the provided function. The method returns a new
        LabeledNeighborPairs object with the specified attribute renamed.

        Args:
            field (str): The attribute to rename.
            func (Callable[[Any], Any]): The function to use for renaming.

        Returns:
            LabeledNeighborPairs: A new LabeledNeighborPairs object with the specified attribute renamed.

        Example:
            ``` py
            >>> np = LabeledNeighborPairs(...)
            >>> np.rename('names', lambda x: x + '_new')
            ```
        """
        return self._apply_func(field, func)

    def intersection(self, other: Self) -> Self:
        """Return the intersection of two LabeledNeighborPairs objects.

        The method calculates the intersection of the pairs from `self` and `other`, and then returns a new
        LabeledNeighborPairs object that contains the intersecting pairs along with their corresponding distances.

        Args:
            other: The other LabeledNeighborPairs object.

        Returns:
            intersected_pairs: A LabeledNeighborPairs object containing the pairs and their corresponding \
                distances that are common between `self` and `other`.

        Example:
            ``` py
            >>> np1 = LabeledNeighborPairs(...)
            >>> np2 = LabeledNeighborPairs(...)
            # 2 ways to get the intersection
            >>> np_intersected = np1 & np2
            >>> np_intersected = np1.intersection(np2)
            ```
        """
        pairs, other_pairs = encode_labels(self.pairs, other.pairs)
        mask = au.intersection(pairs, other_pairs)
        return self.__class__(self.pairs[mask])

    def union(self, other: Self) -> Self:
        """Return the union of two LabeledNeighborPairs objects.

        The method finds the union of the pairs from `self` and `other`. It also ensures that the distances
        in the resulting object correspond to the union pairs.

        Args:
            other: The other LabeledNeighborPairs object to be unified with.

        Returns:
            pairs: A LabeledNeighborPairs object containing the union of the pairs from `self` and `other`, \
                and with corresponding distances.

        Example:
            ``` py
            >>> np1 = LabeledNeighborPairs(...)
            >>> np2 = LabeledNeighborPairs(...)
            # 2 ways to get the union
            >>> np_union = np1 + np2
            >>> np_union = np1.union(np2)
            ```
        """
        pairs, other_pairs = encode_labels(self.pairs, other.pairs)
        mask_a, mask_b = au.union_masks(pairs, other_pairs)

        merged_pairs = np.concatenate((self.pairs[mask_a], other.pairs[mask_b]), axis=0)
        return self.__class__(merged_pairs)

    def difference(self, other: Self) -> Self:
        """Return the difference between two LabeledNeighborPairs objects.

        The method calculates the difference between the pairs from `self` and `other`, then returns a new
        LabeledNeighborPairs object that contains the pairs from `self` that are not in `other`, along with their
        corresponding distances.

        Args:
            other: The other LabeledNeighborPairs object.

        Returns:
            difference_pairs: A LabeledNeighborPairs object containing the pairs and their corresponding \
                distances from `self` that are not in `other`.

        Example:
            ``` py
            >>> np1 = LabeledNeighborPairs(...)
            >>> np2 = LabeledNeighborPairs(...)
            # 2 ways to get the difference
            >>> np_diff = np1 - np2
            >>> np_diff = np1.difference(np2)
            ```

        """
        pairs, other_pairs = encode_labels(self.pairs, other.pairs)
        mask = au.difference(pairs, other_pairs)

        return self.__class__(self.pairs[mask])

    def symmetric_difference(self, other: Self) -> Self:
        """Return the symmetric difference of two LabeledNeighborPairs objects.

        This method creates a new `LabeledNeighborPairs` object that contains pairs and distances
        that are unique to `self` or `other`, but not both.

        Args:
            other: The other LabeledNeighborPairs object.

        Returns:
            A LabeledNeighborPairs object containing the symmetric difference of the two \
                LabeledNeighborPairs objects.

        Example:
            ``` py
            >>> np1 = LabeledNeighborPairs(...)
            >>> np2 = LabeledNeighborPairs(...)
            # 2 ways to get the symmetric difference
            >>> np_sym_diff = np1 | np2
            >>> np_sym_diff = np1.symmetric_difference(np2)
            ```
        """
        pairs, other_pairs = encode_labels(self.pairs, other.pairs)
        mask_a, mask_b = au.symmetric_difference(pairs, other_pairs)

        merged_pairs = np.concatenate((self.pairs[mask_a], other.pairs[mask_b]), axis=0)

        return self.__class__(merged_pairs)

    def isdisjoint(self, other: Self) -> bool:
        """Check if the intersection of two LabeledNeighborPairs objects is null.

        This method checks whether the intersection of the two LabeledNeighborPairs objects is null,
        thus determining if the two objects are disjoint.

        Args:
            other (LabeledNeighborPairs): The other LabeledNeighborPairs object.

        Returns:
            bool: True if the two LabeledNeighborPairs objects are disjoint \
                    (i.e., have no common pairs), and False otherwise.

        Example:
            ``` py
            >>> np1 = LabeledNeighborPairs(...)
            >>> np2 = LabeledNeighborPairs(...)
            >>> np1.isdisjoint(np2)
            ```
            True
        """
        pairs, other_pairs = encode_labels(self.pairs, other.pairs)
        return au.isdisjoint(pairs, other_pairs)

    def issubset(self, other: Self) -> bool:
        """Check if all elements (pairs) of a LabeledNeighborPairs object are found in another LabeledNeighborPairs.

        This method checks whether every pair of atoms from the current LabeledNeighborPairs object
        is also present in the other LabeledNeighborPairs object, \
            thus determining if this object is a subset of the 'other'.

        Args:
            other (LabeledNeighborPairs): The other LabeledNeighborPairs object.

        Returns:
            bool: True if every pair in the current LabeledNeighborPairs object is found
                    in the other LabeledNeighborPairs object, and False otherwise.

        Example:
            ``` py
            >>> np1 = LabeledNeighborPairs(...)
            >>> np2 = LabeledNeighborPairs(...)
            >>> np1.issubset(np2)
            ```
        """
        pairs, other_pairs = encode_labels(self.pairs, other.pairs)
        return au.issubset(pairs, other_pairs)

    def issuperset(self, other: Self) -> bool:
        """Determine if all pairs from another LabeledNeighborPairs object are found in this object.

        This method checks whether every pair of atoms from the 'other' LabeledNeighborPairs object
        is also present in this object, thus determining if this object is a superset of the 'other'.

        Args:
            other (LabeledNeighborPairs): The other LabeledNeighborPairs object.

        Returns:
            bool: True if all pairs from 'other' are found in this object, False otherwise.

        Example:
            ``` py
            >>> np1 = LabeledNeighborPairs(...)
            >>> np2 = LabeledNeighborPairs(...)
            >>> np1.issuperset(np2)
            ```
            True
        """
        pairs, other_pairs = encode_labels(self.pairs, other.pairs)
        return au.issuperset(pairs, other_pairs)

    def isequal(self, other: Type[Self] | object) -> bool:
        """Check if this LabeledNeighborPairs object is equal to another.

        Two LabeledNeighborPairs objects are considered equal if they contain exactly the same pairs.

        Args:
            other (LabeledNeighborPairs): The other LabeledNeighborPairs object.

        Returns:
            bool: True if the two LabeledNeighborPairs objects contain the same pairs, False otherwise.

        Example:
            ``` py
            >>> np1 = LabeledNeighborPairs(...)
            >>> np2 = LabeledNeighborPairs(...)
            >>> np1.isequal(np2)
            True
            ```
        """
        if not isinstance(other, LabeledNeighborPairs):
            warnings.warn(
                f"Comparing {self.__class__.__name__} to {other.__class__.__name__} is not supported.",
                UserWarning,
                stacklevel=2,
            )
            return False
        pairs, other_pairs = encode_labels(self.pairs, other.pairs)
        return au.isequal(pairs, other_pairs)

    def isunique(self) -> bool:
        """Check if all pairs in this LabeledNeighborPairs object are unique.

        This method checks if all pairs in this LabeledNeighborPairs object are unique, i.e., \
            there are no duplicate pairs.

        Returns:
            bool: True if all pairs in this object are unique, False otherwise.

        Example:
            ``` py
            >>> np = LabeledNeighborPairs(...)
            >>> np.isunique()
            False
            ```
        """
        indices, _ = pd.factorize(self.pairs.ravel())  # type: ignore
        pairs = indices.reshape(self.pairs.shape)
        return au.isunique(pairs)

    def is_strict_subset(self, other: Self) -> bool:
        """Check if all pairs of this LabeledNeighborPairs object are in another, and the two sets are not equal.

        A strict subset has all pairs in the 'other' object but the two sets are not identical.

        Args:
            other (LabeledNeighborPairs): The other LabeledNeighborPairs object.

        Returns:
            bool: True if this object is a strict subset of 'other', False otherwise.

        Example:
            ``` py
            >>> np1 = LabeledNeighborPairs(...)
            >>> np2 = LabeledNeighborPairs(...)
            >>> np1.is_strict_subset(np2)
            True
            ```
        """
        pairs, other_pairs = encode_labels(self.pairs, other.pairs)
        return au.is_strict_subset(pairs, other_pairs)

    def is_strict_superset(self, other: Self) -> bool:
        """Check if all pairs of another LabeledNeighborPairs object are in this one, and the two sets are not equal.

        A strict superset has all pairs from the 'other' object but the two sets are not identical.

        Args:
            other (LabeledNeighborPairs): The other LabeledNeighborPairs object.

        Returns:
            bool: True if this object is a strict superset of 'other', False otherwise.

        Example:
            ``` py
            >>> np1 = LabeledNeighborPairs(...)
            >>> np2 = LabeledNeighborPairs(...)
            >>> np1.is_strict_superset(np2)
            True
            ```
        """
        pairs, other_pairs = encode_labels(self.pairs, other.pairs)
        return au.is_strict_superset(pairs, other_pairs)

    @classmethod
    def create_new(cls, pairs: NDArray[np.void]) -> Self:
        """Return a new LabeledNeighborPairs object that is a copy of the current object,
        but with specified pairs.

        Args:
            pairs (NDArray[np.str_]): The labeled pairs to be set.

        Returns:
            A new LabeledNeighborPairs object with the provided pairs.
        """
        cls_instance = cls.__new__(cls)
        cls_instance._pairs = pairs  # noqa: SLF001

        return cls_instance

    def to_frame(self) -> pd.DataFrame:
        """Convert the LabeledNeighborPairs object to a pandas DataFrame.

        Returns:
            A pandas DataFrame containing the atom pairs.
        """
        return pd.DataFrame(self.pairs, columns=["atom1", "atom2"])

    def to_dict(self) -> dict[str, Any]:
        """Convert the LabeledNeighborPairs object to a dictionary.

        Returns:
            A dictionary representation of the LabeledNeighborPairs object.
        """
        return self.to_frame().to_dict(orient="list")  # type: ignore

    @property
    def partner1(self) -> AtomGroupType:
        """Get the first partner of the pairs of indices of atoms that are neighbors.

        Returns:
            The first partner of the atom pairs.
        """
        raise NotImplementedError("This method is not implemented for LabeledNeighborPairs.")

    @property
    def partner2(self) -> AtomGroupType:
        """Get the second partner of the pairs of indices of atoms that are neighbors.

        Returns:
            The second partner of the atom pairs.
        """
        raise NotImplementedError("This method is not implemented for LabeledNeighborPairs.")

    @property
    def pairs(self) -> NDArray[np.void]:
        """Get the pairs of atoms that are neighbors.

        Returns:
            An array containing the pairs of indices of neighboring atoms.
        """
        return self._pairs

    @property
    def labels(self) -> NDArray[np.void]:
        """Get the labels of the pairs of atoms that are neighbors.

        Returns:
            An array containing the labels of the pairs of neighboring atoms.
        """
        return self._pairs

    @property
    def distances(self) -> NDArray[np.float32]:
        """Get the distances between the pairs of indices of atoms that are neighbors.

        Returns:
            An array containing the distances between the pairs of indices of neighboring atoms.
        """
        raise NotImplementedError("This method is not implemented for LabeledNeighborPairs.")

    @property
    def indices(self) -> NDArray[np.int32]:
        """Get the indices of the atoms that are neighbors.

        Returns:
            An array containing the unique indices of the neighboring atoms.
        """
        raise NotImplementedError("This method is not implemented for LabeledNeighborPairs.")

    def __getitem__(self, item: int | slice | NDArray[np.int32]) -> Self:
        """Retrieve the neighbor pairs at the specified index or indices.

        This method allows accessing the neighbor pairs similar to elements in a list.
        For an integer input, it returns a LabeledNeighborPairs object containing a single pair,\
        while for a slice or an array of indices, it returns a LabeledNeighborPairs object with the corresponding pairs.

        Args:
            item (int, slice, or ndarray): An integer, slice, or array of integers indicating
                                            the index/indices of the pair(s).

        Returns:
            LabeledNeighborPairs: A new LabeledNeighborPairs object containing the specified pair(s) and
                            their corresponding distance(s).
        """
        if isinstance(item, int):
            return self.create_new(
                self.pairs[item],
            )

        return self.create_new(self.pairs[item])

    def __contains__(self, other: Self) -> bool:
        """Check whether all pairs in the 'other' LabeledNeighborPairs object are found in this one.

        This method allows using the Python built-in `in` keyword to check for the presence of pairs.
        This method add support for the `in` operator.

        Args:
            other (LabeledNeighborPairs): The LabeledNeighborPairs object to check.

        Returns:
            bool: True if all pairs from the 'other' LabeledNeighborPairs object are found in this one; False otherwise.

        Raises:
            NotImplemented: If the 'other' object is not an instance of the LabeledNeighborPairs class.
        """
        if other.__class__ != self.__class__:
            return NotImplemented

        return au.issubset(other.pairs, self.pairs)

    def __add__(self, other: Self) -> Self:
        """Combine this LabeledNeighborPairs object with another one.

        The resulting LabeledNeighborPairs object is the union of the two sets, containing all unique pairs from both.
        This method adds support for the `+` operator.

        Args:
            other (LabeledNeighborPairs): Another LabeledNeighborPairs object.

        Returns:
            LabeledNeighborPairs: A new LabeledNeighborPairs object that is the union of this one and the 'other'.

        Raises:
            NotImplemented: If the 'other' object is not an instance of the LabeledNeighborPairs class.
        """
        if other.__class__ != self.__class__:
            return NotImplemented
        return self.union(other)

    def __sub__(self, other: Self) -> Self:
        """Get the pairs in this LabeledNeighborPairs object that are not in the 'other'.

        The resulting LabeledNeighborPairs object is the difference of the two sets,
        containing pairs present in this object but not in the 'other'.
        This method adds support for the `-` operator.

        Args:
            other (LabeledNeighborPairs): Another LabeledNeighborPairs object.

        Returns:
            LabeledNeighborPairs: A new LabeledNeighborPairs object that is the difference of this one and the 'other'.

        Raises:
            NotImplemented: If the 'other' object is not an instance of the LabeledNeighborPairs class.
        """
        if other.__class__ != self.__class__:
            return NotImplemented
        return self.difference(other)

    def __or__(self, other: Self) -> Self:
        """Get the pairs that are in either this LabeledNeighborPairs object or the 'other', but not in both.

        The resulting LabeledNeighborPairs object is the symmetric difference of the two sets,
        containing pairs present in either this object or the 'other', but not in both.
        This method adds support for the `|` operator.

        Args:
            other (LabeledNeighborPairs): Another LabeledNeighborPairs object.

        Returns:
            LabeledNeighborPairs: A new LabeledNeighborPairs object that is the \
                symmetric difference of this one and the 'other'.

        Raises:
            NotImplemented: If the 'other' object is not an instance of the LabeledNeighborPairs class.
        """
        if other.__class__ != self.__class__:
            return NotImplemented
        return self.symmetric_difference(other)

    def __eq__(self, other: object) -> bool:
        """Check if this LabeledNeighborPairs object is equal to the 'other'.

        Equality is based on the pairs and their distances.
        This method adds support for the `==` operator.

        Args:
            other (Any): Another object.

        Returns:
            bool: True if 'other' is an identical LabeledNeighborPairs object; False otherwise.

        Raises:
            NotImplemented: If the 'other' object is not an instance of the LabeledNeighborPairs class.
        """
        if other.__class__ != self.__class__:
            return NotImplemented
        return self.isequal(other)

    def __and__(self, other: Self) -> Self:
        """Get the pairs that are common to both this LabeledNeighborPairs object and the 'other'.

        The resulting LabeledNeighborPairs object is the intersection of the two sets,
        containing pairs present in both this object and the 'other'.
        This method adds support for the `&` operator.

        Args:
            other (LabeledNeighborPairs): Another LabeledNeighborPairs object.

        Returns:
            LabeledNeighborPairs: A new LabeledNeighborPairs object that is the \
                intersection of this one and the 'other'.

        Raises:
            NotImplemented: If the 'other' object is not an instance of the LabeledNeighborPairs class.
        """
        if other.__class__ != self.__class__:
            return NotImplemented
        return self.intersection(other)

    def __xor__(self, other: Self) -> Self:
        """Get the pairs that are in either this LabeledNeighborPairs object or the 'other', but not in both.

        This method behaves similarly to the `__or__` method.
        This method adds support for the `^` operator.

        Args:
            other (LabeledNeighborPairs): Another LabeledNeighborPairs object.

        Returns:
            LabeledNeighborPairs: A new LabeledNeighborPairs object that is the \
                symmetric difference of this one and the 'other'.

        Raises:
            NotImplemented: If the 'other' object is not an instance of the LabeledNeighborPairs class.
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

    def __ne__(self, other: object) -> bool:
        if other.__class__ != self.__class__:
            return NotImplemented
        return not self.isequal(other)

    def __len__(self) -> int:
        """Get the number of pairs in this LabeledNeighborPairs object."""
        return self.pairs.shape[0]

    def __str__(self) -> str:
        return f"<Lahuta LabeledNeighborPairs class containing {self.pairs.shape[0]} pairs>"

    def __repr__(self) -> str:
        return self.__str__()
