"""
Placeholder
"""
from __future__ import annotations

from typing import Tuple, TypeVar

import numpy as np
import numpy.typing as npt
from numpy.typing import NDArray

# from typing_extensions import TypeVarTuple, Unpack

_DType = TypeVar("_DType", np.float32, np.int32)
NDArrayDType = NDArray[_DType]
T1 = TypeVar("T1", bound=npt.NBitBase)
T2 = TypeVar("T2", bound=npt.NBitBase)
NDArrayInt = npt.NDArray[np.int32]

# _NewShape = TypeVarTuple("_NewShape")
# # Represents an NDArray reshaped with np.newaxis
# # ReshapedArray = NewType('ReshapedArray', Tuple[NDArray[np.float32], _NewShape])
# ReshapedArray = NewType("ReshapedArray", Tuple[NDArray[np.float32], Unpack[_NewShape]])


# pylint: disable=W1114
def calculate_angle(
    vertex: NDArray[np.float32],
    point1: NDArray[np.float32],
    point2: NDArray[np.float32],
    degrees: bool = False,
) -> NDArray[np.float32]:
    """Calculate the angle between three points.

    This function calculates the angle at `vertex` created by points `point1` and `point2`.

    Args:
        vertex (NDArray[np.float32]): The coordinates of the vertex point, where the angle is being measured.
        point1 (NDArray[np.float32]): The coordinates of the first point.
        point2 (NDArray[np.float32]): The coordinates of the second point.
        degrees (bool, optional): If True, the angle is returned in degrees. Otherwise, it's returned in radians.
            Default is False (radians).

    Returns:
        angle (float): The calculated angle in radians (default) or degrees, depending on the value of the `degrees` argument.

    """
    # Vector from vertex to point1 & point2
    vector1: NDArray[np.float32] = point1 - vertex
    vector2: NDArray[np.float32] = point2 - vertex

    # Normalize vector1
    magnitude_vector1: NDArray[np.float32] = np.linalg.norm(vector1, axis=-1)
    normalized_vector1: NDArray[np.float32] = (
        vector1 / magnitude_vector1[:, :, np.newaxis]
    )

    # Normalize vector2
    magnitude_vector2: NDArray[np.float32] = np.linalg.norm(vector2, axis=-1)
    normalized_vector2: NDArray[np.float32] = (
        vector2 / magnitude_vector2[:, :, np.newaxis]
    )

    # Dot product of normalized vectors
    dot_product: NDArray[np.float32] = np.sum(
        normalized_vector1 * normalized_vector2, axis=-1
    )

    angle_rad: NDArray[np.float32] = np.arccos(dot_product)

    if degrees:
        return np.degrees(angle_rad)

    return angle_rad


def array_distance(
    arr1: NDArray[np.float32], arr2: NDArray[np.float32]
) -> NDArray[np.float32]:
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
    arr_reshape = arr1[:, np.newaxis, :]
    distance_array: NDArray[np.float32] = np.linalg.norm(arr_reshape - arr2, axis=-1)

    return distance_array


def matching_indices(
    arr1: NDArray[np.int32], arr2: NDArray[np.int32]
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

    arr1_reshape = arr1[:, None]
    idx: NDArray[np.bool_] = (arr1_reshape == arr2).all(-1).any(1)
    return idx


def optimized_matching_pairs(
    arr1: NDArray[np.int32], arr2: NDArray[np.int32]
) -> NDArray[np.int32]:
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
    set2 = set(map(tuple, arr2))

    # Find common elements
    common = np.array(list(set1 & set2))

    return common


def np_optimized_matching_pairs(
    arr1: NDArray[np.int32], arr2: NDArray[np.int32]
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
    arr1: NDArray[np.int32], arr2: NDArray[np.int32]
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
    idx: NDArray[np.bool_] = (arr1[:, None] != arr2).any(-1).all(1)

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

    concat = np.concatenate((arr1, arr2), axis=0)

    unique_indices = np.unique(concat, axis=0, return_index=True)[1]
    sorted_indices = np.sort(unique_indices)

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
