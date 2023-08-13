"""
module: `lahuta.utils.array_utils.py`

This module contains a suite of utility functions for operations on 2D numpy arrays, particularly arrays 
representing atom pairs. These operations include finding shared pairs, determining non-matching indices,
and various set-like operations including intersection, difference, symmetric difference, union, and 
subset/superset/equality checks. It also provides functions to check for uniqueness 
and strict subset/superset relations.

Functions:
    ```
    find_shared_pairs(arr1, arr2): Finds shared pairs between two arrays.
    non_matching_indices(arr1, arr2): Finds non-matching indices between two arrays.
    intersection(arr1, arr2): Computes the intersection of two arrays.
    difference(arr1, arr2): Computes the difference of two arrays.
    symmetric_difference(arr1, arr2): Computes the symmetric difference of two arrays.
    union(arr1, arr2): Computes the union of two arrays.
    isdisjoint(arr1, arr2): Checks if two arrays are disjoint.
    issubset(arr1, arr2): Checks if the first array is a subset of the second array.
    issuperset(arr1, arr2): Checks if the first array is a superset of the second array.
    isequal(arr1, arr2): Checks if two arrays are equal.
    isunique(arr): Checks if an array has unique rows.
    is_strict_subset(arr1, arr2): Checks if the first array is a strict subset of the second array.
    is_strict_superset(arr1, arr2): Checks if the first array is a strict superset of the second array.
    ```

Each function operates on numpy arrays (with a shape of (n,2) for most functions,
representing pairs of atom indices) and returns either a new array resulting from the 
operation or a boolean value representing the relationship between arrays.

Notes:
    This module is intended for use with arrays of atom pair indices. However, most of these 
    functions would be applicable to other data as long as the input is 2D numpy arrays.

    Functions like `intersection`, `difference`, `union` etc. perform set operations considering each row of the input 
    arrays as an element of the set. This makes these functions particularly useful for operations on collections of 
    atom pairs, where each pair is represented by a row in the array.
"""

from typing import Tuple, TypeVar

import numpy as np
import numpy.typing as npt
from numpy.typing import NDArray

# from typing_extensions import TypeVarTuple, Unpack

_DType = TypeVar("_DType", np.float32, np.int32, np.void)
NDArrayDType = NDArray[_DType]
NDArrayInt = npt.NDArray[np.int32]

__all__ = [
    "find_shared_pairs",
    "non_matching_indices",
    "intersection",
    "difference",
    "symmetric_difference",
    "union",
    "isdisjoint",
    "issubset",
    "issuperset",
    "isequal",
    "isunique",
    "is_strict_subset",
    "is_strict_superset",
]


def sorting_indices(arr: NDArray[np.int32]) -> NDArray[np.int32]:
    """
    Sort pairs in an array such that the smaller index is in the first column.

    Args:
        arr (NDArray[np.int32]): The array to sort.

    Returns:
        NDArray[np.int32]: The sorted array.
    """
    # Ensure the smaller index is in the first column
    arr = np.sort(arr, axis=1)

    # Use lexsort to get sorted indices from large to small,
    # then use it to index into the sorted array
    indices: NDArray[np.int32] = np.lexsort((arr[:, 1], arr[:, 0]))  # type: ignore

    return indices


def sort_pairs(arr: NDArray[np.int32]) -> NDArray[np.int32]:
    """
    Sort pairs in an array such that the smaller index is in the first column.

    Args:
        arr (NDArray[np.int32]): The array to sort.

    Returns:
        NDArray[np.int32]: The sorted array.
    """
    indices = sorting_indices(arr)

    return arr[indices]


def find_shared_pairs(arr1: NDArray[np.int32], arr2: NDArray[np.int32]) -> NDArray[np.bool_]:
    """Find shared elements between two 2D numpy arrays.

    This function takes two 2D arrays where each row represents a pair of atom indices and returns a 1D boolean
    array representing whether each pair in `arr1` also appears in `arr2`.

    Args:
        arr1 (NDArray[np.int32]): A 2D array of shape (n_pairs1, 2) where each row represents a pair of atom indices.
        arr2 (NDArray[np.int32]): A 2D array of shape (n_pairs2, 2) where each row represents a pair of atom indices.

    Returns:
        NDArray[np.bool_]: A 1D boolean array of shape (n_pairs1,) where each element represents
        whether the corresponding pair in `arr1` appears in `arr2`.

    ???+ example "Example"
        ``` py
        arr1 = np.array([[1, 2], [3, 4], [5, 6]])
        arr2 = np.array([[3, 4], [7, 8], [1, 2]])
        find_shared_pairs(arr1, arr2)
        array([ True,  True, False])
        ```
    """
    arr1_complex = arr1[:, 0] + 1j * arr1[:, 1]
    arr2_complex = arr2[:, 0] + 1j * arr2[:, 1]
    return np.isin(arr1_complex, arr2_complex)


def non_matching_indices(arr1: NDArray[np.int32], arr2: NDArray[np.int32]) -> NDArray[np.bool_]:
    """
    Find the indices of non-matching elements between two 2D numpy arrays.

    This function takes two 2D arrays where each row represents a pair of atom indices and returns a 1D boolean
    array representing whether each pair in `arr1` does not appear in `arr2`.

    Args:
        arr1 (NDArray[np.int32]): A 2D array of shape (n_pairs1, 2) where each row represents a pair of atom indices.

    Returns:
        NDArray[np.bool_]: A 1D boolean array of shape (n_pairs1,) where each element
        represents whether the corresponding

    ???+ example "Example"
        ``` py
        arr1 = np.array([[1, 2], [3, 4], [5, 6]])
        arr2 = np.array([[3, 4], [7, 8], [1, 2]])
        non_matching_indices(arr1, arr2)
        array([False, False,  True])
        ```
    """
    idx: NDArray[np.bool_] = (arr1[:, None] != arr2).any(-1).all(1)

    return idx


# def asvoid(arr: NDArray[_DType]) -> NDArray[np.void]:
#     """
#     Reference: https://stackoverflow.com/questions/16216078/test-for-membership-in-a-2d-numpy-array
#     Based on http://stackoverflow.com/a/16973510/190597 (Jaime, 2013-06)
#     View the array as dtype np.void (bytes). The items along the last axis are
#     viewed as one value. This allows comparisons to be performed on the entire row.
#     """
#     arr = np.ascontiguousarray(arr)
#     if np.issubdtype(arr.dtype, np.floating):
#         # Care needs to be taken here since
#         # np.array([-0.]).view(np.void) != np.array([0.]).view(np.void)
#         # Adding 0. converts -0. to 0.
#         arr += 0.0  # type: ignore
#     return arr.view(np.dtype((np.void, arr.dtype.itemsize * arr.shape[-1])))  # type: ignore


def asvoid(arr: NDArray[_DType]) -> NDArray[np.void]:
    """
    Convert a 2D numpy array into a 1D array of type np.void.

    The function views each row in the input array as a single item of type np.void, which effectively
    converts the 2D array into a 1D array of binary representations of the rows. This can be useful
    for performing operations that are typically only possible with 1D arrays, such as membership tests.

    Args:
        arr: A 2D numpy array of shape (n, 2), where each row is a pair of atom indices.

    Returns:
        arr_void: A 1D numpy array of dtype np.void, where each element is the binary representation
                  of a pair of atom indices from the input array.

    ???+ example "Example"
        ``` py
        arr = np.array([[1, 2], [3, 4]])
        print(asvoid(arr))
        [(1, 2), (3, 4)]
        ```
    """
    arr = np.ascontiguousarray(arr)
    if np.issubdtype(arr.dtype, np.floating):
        arr += 0.0  # type: ignore
    return arr.view(np.dtype((np.void, arr.dtype.itemsize * arr.shape[-1])))  # type: ignore


def intersection(arr1: NDArray[_DType], arr2: NDArray[_DType], assume_unique: bool = False) -> NDArray[np.bool_]:
    """
    Calculate the intersection of two 2D arrays and return a boolean array that indicates the membership
    of each pair of atom indices in `arr1` in `arr2`.

    The function first converts `arr1` and `arr2` into 1D arrays of type np.void using the `asvoid` function,
    then uses `np.in1d` to test the membership of each element in the resulting 1D array from `arr1` in that
    from `arr2`.

    Args:
        arr1: A 2D numpy array of shape (n, 2) where each row is a pair of atom indices.
        arr2: A 2D numpy array of shape (m, 2) where each row is a pair of atom indices.
        assume_unique: If True, the input arrays are both assumed to be unique, which can speed up the
                       calculation. Default is False.

    Returns:
        mask: A 1D boolean numpy array of length n, where each element indicates whether the corresponding
              pair of atom indices in `arr1` is also in `arr2`.

    ???+ example "Example"
        ``` py
        arr1 = np.array([[1, 2], [3, 4]])
        arr2 = np.array([[1, 2], [5, 6]])
        print(intersection(arr1, arr2))
        [True False]
        ```
    """
    arr1_void = asvoid(arr1)
    arr2_void = asvoid(arr2)
    return np.in1d(arr1_void, arr2_void, assume_unique)


def difference(arr1: NDArray[_DType], arr2: NDArray[_DType], assume_unique: bool = False) -> NDArray[np.bool_]:
    """
    Calculate the set difference between two 2D arrays, represented as a boolean array.

    The function converts the input arrays into 1D arrays of type np.void, then uses `np.in1d` to test
    the membership of each element in `arr1` in `arr2`. The resulting boolean array is then inverted to
    represent the elements in `arr1` that are not in `arr2`.

    Args:
        arr1: A 2D numpy array of shape (n, 2) where each row is a pair of atom indices.
        arr2: A 2D numpy array of shape (m, 2) where each row is a pair of atom indices.
        assume_unique: If True, the input arrays are both assumed to be unique, which can speed up the
                       calculation. Default is False.

    Returns:
        mask: A 1D boolean numpy array of length n, where each element indicates whether the corresponding
              pair of atom indices in `arr1` is not in `arr2`.

    ???+ example "Example"
        ``` py
        arr1 = np.array([[1, 2], [3, 4]])
        arr2 = np.array([[1, 2], [5, 6]])
        print(difference(arr1, arr2))
        [False True]
        ```
    """

    arr1_void = asvoid(arr1)
    arr2_void = asvoid(arr2)

    return np.in1d(arr1_void, arr2_void, assume_unique, invert=True)  # type: ignore


def symmetric_difference(
    arr1: NDArray[_DType], arr2: NDArray[_DType], assume_unique: bool = False
) -> Tuple[NDArray[np.bool_], NDArray[np.bool_]]:
    """Calculate the symmetric difference of two arrays.

    This function returns the elements that are in `arr1` but not in `arr2` and vice versa.

    Args:
        arr1: An array of shape (n, 2) where each row is a pair of atom indices.
        arr2: An array of shape (m, 2) where each row is a pair of atom indices.
        assume_unique: If True, the input arrays are both assumed to be unique,
                       which can speed up the calculation. Default is False.

    Returns:
        mask_a: A boolean array that can be used to index `arr1` to get the elements unique to `arr1`.
        mask_b: A boolean array that can be used to index `arr2` to get the elements unique to `arr2`.
    """
    # pylint: disable=arguments-out-of-order
    mask_a = difference(arr1, arr2, assume_unique)
    mask_b = difference(arr2, arr1, assume_unique)

    return mask_a, mask_b


def union(arr1: NDArray[_DType], arr2: NDArray[_DType]) -> Tuple[NDArray[_DType], NDArray[np.int32]]:
    """
    Calculate the union of two arrays and return the unique pairs along with their indices.

    The function finds unique pairs from the union of `arr1` and `arr2`, and also returns the
    indices in the concatenated array which correspond to these unique pairs. The indices can
    be used to select corresponding elements in another array that has the same order as `arr1`
    and `arr2`.

    Args:
        arr1: A 2D numpy array of shape (n, 2) where each row is a pair of atom indices.
        arr2: A 2D numpy array of shape (m, 2) where each row is a pair of atom indices.

    Returns:
        union_arr: A 2D numpy array of shape (k, 2) where each row is a unique pair of atom indices
                   that are in `arr1` or `arr2`.
        indices: A 1D numpy array of the indices in the concatenated array which correspond to the
                 unique pairs.

    ???+ example "Example"
        ``` py
        arr1 = np.array([[1, 2], [3, 4]])
        arr2 = np.array([[1, 2], [5, 6]])
        pairs, indices = union(arr1, arr2)
        print(pairs)
        print(indices)
        [[1, 2], [3, 4], [5, 6]]
        [0 1 3]
        ```
    """
    concat = np.concatenate((arr1, arr2), axis=0)

    _, unique_indices = np.unique(concat, axis=0, return_index=True)
    sorted_indices = np.sort(unique_indices)

    return concat[sorted_indices], sorted_indices


def union_masks(
    arr1: NDArray[_DType], arr2: NDArray[_DType], assume_unique: bool = False
) -> Tuple[NDArray[np.bool_], NDArray[np.bool_]]:
    """Calculate the union of two arrays.

    Args:
        arr1: An array of shape (n, 2) where each row is a pair of atom indices.
        arr2: An array of shape (m, 2) where each row is a pair of atom indices.
        assume_unique: If True, the input arrays are both assumed to be unique,
                       which can speed up the calculation. Default is False.

    Returns:
        mask_a: A boolean array that can be used to index `arr1` to get the elements of the union from `arr1`.
        mask_b: A boolean array that can be used to index `arr2` to get the elements of the union from `arr2`.
    """
    # pylint: disable=arguments-out-of-order
    mask_a_diff, mask_b_diff = symmetric_difference(arr1, arr2, assume_unique)
    mask_a_int = intersection(arr1, arr2, assume_unique)

    mask_a_union = mask_a_diff | mask_a_int
    mask_b_union = mask_b_diff 

    return mask_a_union, mask_b_union


def isdisjoint(arr1: NDArray[_DType], arr2: NDArray[_DType]) -> bool:
    """
    Determines if two arrays have a null intersection.

    This function checks if the two input arrays do not share any common elements. Both arrays are assumed
    to be 2D with each row being a pair of atom indices.

    Args:
        arr1 (np.ndarray): A 2D numpy array of shape (n, 2) where each row is a pair of atom indices.
        arr2 (np.ndarray): A 2D numpy array of shape (m, 2) where each row is a pair of atom indices.

    Returns:
        np.bool_: True if the two arrays have no common elements (i.e., disjoint), False otherwise.

    ???+ example "Example"
        ``` py
        arr1 = np.array([[1, 2], [3, 4]])
        arr2 = np.array([[5, 6], [7, 8]])
        isdisjoint(arr1, arr2)
        True
        ```
    """

    # return arr1[intersection(arr1, arr2)].size == 0
    return not np.any(intersection(arr1, arr2))


def issubset(arr1: NDArray[_DType], arr2: NDArray[_DType]) -> bool:
    """
    Determines if the first array is a subset of the second array.

    This function checks if all elements of the first array are found within the second array. Both arrays are
    assumed to be 2D with each row being a pair of atom indices.

    Args:
        arr1 (np.ndarray): A 2D numpy array of shape (n, 2) where each row is a pair of atom indices.
        arr2 (np.ndarray): A 2D numpy array of shape (m, 2) where each row is a pair of atom indices.

    Returns:
        bool: True if all elements of the first array are in the second array
              (i.e., arr1 is a subset of arr2), False otherwise.

    ???+ example "Example"
        ``` py
        arr1 = np.array([[1, 2], [3, 4]])
        arr2 = np.array([[1, 2], [3, 4], [5, 6]])
        issubset(arr1, arr2)
        True
        ```
    """

    # return bool(np.all(intersection(arr1, arr2)))  # type: ignore
    result = np.sum(intersection(arr1, arr2)) == len(arr1)
    return bool(result)


def issuperset(arr1: NDArray[_DType], arr2: NDArray[_DType]) -> bool:
    """
    Checks if `arr1` is a superset of `arr2`.

    This function takes two 2D numpy arrays, arr1 and arr2, each with shape (n, 2),
    and checks if every pair in arr2 is present in arr1. Each row in the arrays represents
    a pair of atom indices. The function returns True if arr1 is a superset of arr2,
    and False otherwise.

    Args:
        arr1 (np.ndarray): A 2D numpy array of shape (n, 2) representing pairs of atom indices.
        arr2 (np.ndarray): A 2D numpy array of shape (n, 2) representing pairs of atom indices.

    Returns:
        bool: True if `arr1` is a superset of `arr2`, False otherwise.

    ???+ example "Example"
        ``` py
        issuperset(np.array([[1, 2], [2, 3]]), np.array([[1, 2]]))
        True
        issuperset(np.array([[1, 2], [2, 3]]), np.array([[1, 2], [3, 4]]))
        False
        ```
    """
    # pylint: disable=arguments-out-of-order
    return issubset(arr2, arr1)


def isequal(arr1: NDArray[_DType], arr2: NDArray[_DType]) -> bool:
    """
    Checks if two 2D arrays, `arr1` and `arr2`, are equal.

    This function verifies the equality of two 2D numpy arrays, `arr1` and `arr2`.
    Each array has a shape (n, 2), where each row represents a pair of atom indices.
    The arrays are considered equal if they have the same shape and each pair from
    `arr1` is found in `arr2`, and vice versa.

    Args:
        arr1 (np.ndarray): A 2D numpy array of shape (n, 2) representing pairs of atom indices.
        arr2 (np.ndarray): A 2D numpy array of shape (n, 2) representing pairs of atom indices.

    Returns:
        bool: True if `arr1` and `arr2` are equal, False otherwise.

    ???+ example "Example"
        ``` py
        isequal(np.array([[1, 2], [2, 3]]), np.array([[1, 2], [2, 3]]))
        True
        isequal(np.array([[1, 2], [2, 3]]), np.array([[1, 2], [3, 4]]))
        False
        ```
    """
    if arr1.shape != arr2.shape:
        return False

    return issubset(arr1, arr2) and issubset(arr2, arr1)


def isunique(arr: NDArray[_DType]) -> bool:
    """
    Checks if a 2D array, `arr`, contains no duplicate entries.

    The function checks if a 2D numpy array, `arr`, of shape (n, 2), has any duplicate rows.
    Each row represents a pair of atom indices.
    The function returns True if `arr` contains no duplicate pairs, and False otherwise.

    Args:
        arr (np.ndarray): A 2D numpy array of shape (n, 2) representing pairs of atom indices.

    Returns:
        bool: True if `arr` contains no duplicate pairs, False otherwise.

    ???+ example "Example"
        ``` py
        isunique(np.array([[1, 2], [2, 3]]))
        True
        isunique(np.array([[1, 2], [2, 3], [1, 2]]))
        False
        ```
    """

    return arr.shape[0] == np.unique(arr, axis=0).shape[0]  # type: ignore


def is_strict_subset(arr1: NDArray[_DType], arr2: NDArray[_DType]) -> bool:
    """
    Check if the first array is a strict subset of the second array.

    The function determines whether every element of the first array is in the second array and the two arrays
    are not identical. The arrays should be of shape (n, 2), where each row is a pair of atom indices.

    Args:
        arr1 (np.ndarray): An array of shape (n, 2) where each row is a pair of atom indices.
        arr2 (np.ndarray): An array of shape (n, 2) where each row is a pair of atom indices.

    Returns:
        bool: True if the first array is a strict subset of the second array, False otherwise.

    ???+ example "Example"
        ``` py
        arr1 = np.array([[1, 2], [3, 4]])
        arr2 = np.array([[1, 2], [3, 4], [5, 6]])
        is_strict_subset(arr1, arr2)
        True
        ```
    """

    return issubset(arr1, arr2) and not isequal(arr1, arr2)


def is_strict_superset(arr1: NDArray[_DType], arr2: NDArray[_DType]) -> bool:
    """
    Check if the first array is a strict superset of the second array.

    The function determines whether every element of the second array is in the first array and the two arrays
    are not identical. The arrays should be of shape (n, 2), where each row is a pair of atom indices.

    Args:
        arr1 (np.ndarray): An array of shape (n, 2) where each row is a pair of atom indices.
        arr2 (np.ndarray): An array of shape (n, 2) where each row is a pair of atom indices.

    Returns:
        bool: True if the first array is a strict superset of the second array, False otherwise.

    ???+ example "Example"
        ``` py
        arr1 = np.array([[1, 2], [3, 4], [5, 6]])
        arr2 = np.array([[1, 2], [3, 4]])
        is_strict_superset(arr1, arr2)
        True
        ```
    """

    return issuperset(arr1, arr2) and not isequal(arr1, arr2)
