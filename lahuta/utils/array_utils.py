"""
Placeholder
"""
from typing import Set, Tuple, TypeVar, Union

import numpy as np
from numpy.typing import NDArray

_DType = TypeVar("_DType", np.float_, np.int_)


# pylint: disable=W1114
def calculate_angle(
    point_a: NDArray[np.float_],
    point_b: NDArray[np.float_],
    point_c: NDArray[np.float_],
    degrees: bool = False,
) -> float:
    """Calculate the angle between three points.

    Parameters
    ----------
    point_a : np.ndarray
        The first point.
    point_b : np.ndarray
        The second point.
    point_c : np.ndarray
        The third point.

    Returns
    -------
    angle : float
        The angle in radians.
    """
    v1 = point_a - point_b
    v2 = point_c - point_b

    v1_mag = np.linalg.norm(v1, axis=-1)  # type: ignore
    v1_norm = v1 / v1_mag[:, :, np.newaxis]

    v2_mag = np.linalg.norm(v2, axis=-1)  # type: ignore
    v2_norm = v2 / v2_mag[:, :, np.newaxis]

    res = np.sum(v1_norm * v2_norm, axis=-1)  # type: ignore

    if degrees:
        return np.degrees(np.arccos(res))

    return np.arccos(res)


def array_distance(
    arr1: NDArray[np.float_], arr2: NDArray[np.float_]
) -> NDArray[np.float_]:
    """Takes the difference between the two arrays and calculates the norm of the difference.


    Parameters
    ----------
    arr1: np.array
        Shape (n, 3) array

    arr2: np.array
        Shape (m, 6, 3) array

    Returns
    -------
    distance_array: np.array
        Shape (n, 6) array
    """
    distance_array = np.linalg.norm(arr1[:, np.newaxis, :] - arr2, axis=-1)  # type: ignore

    return distance_array


def matching_indices(
    arr1: NDArray[np.int_], arr2: NDArray[np.int_]
) -> NDArray[np.bool_]:
    """Return indices of elements in `arr1` that are in `arr2`.

    Parameters
    ----------
    arr1 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    arr2 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.

    Returns
    -------
    arr : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.

    """

    idx = (arr1[:, None] == arr2).all(-1).any(1)
    return idx


def optimized_matching_pairs(
    arr1: NDArray[np.int_], set2: Union[NDArray[np.int_], Set[Tuple[int, int]]]
) -> NDArray[np.int_]:
    """Return elements in `arr1` that are in `arr2`.

    Parameters
    ----------
    arr1 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    arr2 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.

    Returns
    -------
    arr : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.

    """
    # Convert arrays to sets of tuples
    set1 = set(map(tuple, arr1))
    if not isinstance(set2, set):
        set2 = set(map(tuple, set2))

    # Find common elements
    common = np.array(list(set1 & set2))

    return common


def np_optimized_matching_pairs(
    arr1: NDArray[np.int_], arr2: NDArray[np.int_]
) -> NDArray[np.bool_]:
    """Return elements in `arr1` that are in `arr2`.

    Parameters
    ----------
    arr1 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    arr2 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.

    Returns
    -------
    arr : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    """

    # Convert the pairs in both arrays to complex numbers
    arr1_complex = arr1[:, 0] + 1j * arr1[:, 1]
    arr2_complex = arr2[:, 0] + 1j * arr2[:, 1]

    # Get common complex numbers
    common_complex = np.in1d(arr1_complex, arr2_complex)  # type: ignore

    return common_complex


def non_matching_indices(
    arr1: NDArray[np.int_], arr2: NDArray[np.int_]
) -> NDArray[np.bool_]:
    """Return the elements in `arr1` that are not in `arr2`.

    Parameters
    ----------
    arr1 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    arr2 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.

    Returns
    -------
    arr : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.

    """
    idx = (arr1[:, None] != arr2).any(-1).all(1)

    return idx


# These are all set operations we want to implement in our classes
def asvoid(arr: NDArray[_DType]) -> NDArray[np.void]:
    """
    Reference: https://stackoverflow.com/questions/16216078/test-for-membership-in-a-2d-numpy-array
    Based on http://stackoverflow.com/a/16973510/190597 (Jaime, 2013-06)
    View the array as dtype np.void (bytes). The items along the last axis are
    viewed as one value. This allows comparisons to be performed on the entire row.
    """
    arr = np.ascontiguousarray(arr)
    if np.issubdtype(arr.dtype, np.floating):
        # Care needs to be taken here since
        # np.array([-0.]).view(np.void) != np.array([0.]).view(np.void)
        # Adding 0. converts -0. to 0.
        arr += 0.0  # type: ignore
    return arr.view(np.dtype((np.void, arr.dtype.itemsize * arr.shape[-1])))  # type: ignore


def intersection(
    arr1: NDArray[_DType], arr2: NDArray[_DType], assume_unique: bool = False
) -> NDArray[np.bool_]:
    """Calculate the intersection of two arrays and return the indices of the
    elements in `arr1` that are in `arr2`.

    Parameters
    ----------
    arr1 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    arr2 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    Returns
    -------
    arr : np.ndarray
        An array of shape (n, 2) where each row is a pair of `arr1` atom indices that are in both `arr1` and `arr2`.
    """

    arr1_void = asvoid(arr1)
    arr2_void = asvoid(arr2)
    return np.in1d(arr1_void, arr2_void, assume_unique)  # type: ignore


def difference(
    arr1: NDArray[_DType], arr2: NDArray[_DType], assume_unique: bool = False
) -> NDArray[np.bool_]:
    """Calculate the difference of two arrays and return the indices of the
    elements in `arr1` that are not in `arr2`.

    Parameters
    ----------
    arr1 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    arr2 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    Returns
    -------
    arr : np.ndarray
        An array of shape (n, 2) where each row is a pair of `arr1` atom indices that are in `arr1` but not in `arr2`.
    """

    arr1_void = asvoid(arr1)
    arr2_void = asvoid(arr2)

    return np.in1d(arr1_void, arr2_void, assume_unique, invert=True)  # type: ignore


def symmetric_difference(
    arr1: NDArray[_DType], arr2: NDArray[_DType], assume_unique: bool = False
) -> Tuple[NDArray[np.bool_], NDArray[np.bool_]]:
    """Calculate the symmetric difference of two arrays and return the indices of the
    elements in `arr1` that are not in `arr2` and the indices of the elements in `arr2` that are not in `arr1`.

    Parameters
    ----------
    arr1 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    arr2 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    Returns
    -------
    arr : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices that are in either `arr1` or `arr2` but not both.
    """

    mask_a = difference(arr1, arr2, assume_unique)
    mask_b = difference(arr2, arr1, assume_unique)

    return mask_a, mask_b


def union(arr1: NDArray[_DType], arr2: NDArray[_DType]) -> NDArray[np.int_]:
    """Calculate the union of two arrays and return the indices of the
    elements in `arr1` and `arr2`. Duplicate entries are removed. Neighbors indices are sorted.


    Parameters
    ----------
    arr1 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    arr2 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    Returns
    -------
    arr : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices that are in `arr1` or `arr2`.
    """

    concat = np.concatenate((arr1, arr2), axis=0)  # type: ignore

    unique_indices = np.unique(concat, axis=0, return_index=True)[1]  # type: ignore
    sorted_indices = np.sort(unique_indices)  # type: ignore

    return sorted_indices


def isdisjoint(arr1: NDArray[_DType], arr2: NDArray[_DType]) -> bool:
    """Return True if two arrays have a null intersection.

    Parameters
    ----------
    arr1 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    arr2 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    Returns
    -------
    bool
        True if the arrays have null intersection.
    """

    return arr1[intersection(arr1, arr2)].size == 0


def issubset(arr1: NDArray[_DType], arr2: NDArray[_DType]) -> bool:
    """Return True if the first array is a subset of the second array.

    Parameters
    ----------
    arr1 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    arr2 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    Returns
    -------
    bool
        True if the first array is a subset of the second array.
    """

    return bool(np.all(intersection(arr1, arr2)))  # type: ignore


def issuperset(arr1: NDArray[_DType], arr2: NDArray[_DType]) -> bool:
    """Return True if the first array is a superset of the second array.

    Parameters
    ----------
    arr1 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    arr2 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    Returns
    -------
    bool
        True if the first array is a superset of the second array.
    """

    return issubset(arr2, arr1)


def isequal(arr1: NDArray[_DType], arr2: NDArray[_DType]) -> bool:
    """Return True if the two arrays are equal.

    Parameters
    ----------
    arr1 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    arr2 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    Returns
    -------
    bool
        True if the two arrays are equal.
    """

    return issubset(arr1, arr2) and issubset(arr2, arr1)


def isunique(arr: NDArray[_DType]) -> bool:
    """Return True if the array contains no duplicate entries.

    Parameters
    ----------
    arr : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    Returns
    -------
    bool
        True if the array contains no duplicate entries.
    """

    return arr.shape[0] == np.unique(arr, axis=0).shape[0]  # type: ignore


def is_strict_subset(arr1: NDArray[_DType], arr2: NDArray[_DType]) -> bool:
    """Return True if the first array is a strict subset of the second array,
    but not identical.

    Parameters
    ----------
    arr1 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    arr2 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    Returns
    -------
    bool
        True if the first array is a strict subset of the second array.
    """

    return issubset(arr1, arr2) and not isequal(arr1, arr2)


def is_strict_superset(arr1: NDArray[_DType], arr2: NDArray[_DType]) -> bool:
    """Return True if the first array is a strict superset of the second array,
    but not identical.

    Parameters
    ----------
    arr1 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    arr2 : np.ndarray
        An array of shape (n, 2) where each row is a pair of atom indices.
    Returns
    -------

    bool
        True if the first array is a strict superset of the second array.
    """

    return issuperset(arr1, arr2) and not isequal(arr1, arr2)
