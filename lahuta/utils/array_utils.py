"""
Placeholder
"""

from functools import partial, update_wrapper

import numpy as np


class ArraySetOps:
    def __init__(self, arr):
        self.arr = arr

        self.intersection = partial(intersection, self.arr)
        self.union = partial(union, self.arr)
        self.difference = partial(difference, self.arr)
        self.symmetric_difference = partial(symmetric_difference, self.arr)
        self.isdisjoint = partial(isdisjoint, self.arr)
        self.issubset = partial(issubset, self.arr)
        self.issuperset = partial(issuperset, self.arr)
        self.isequal = partial(isequal, self.arr)
        self.isunique = partial(isunique, self.arr)
        self.is_strict_subset = partial(is_strict_subset, self.arr)
        self.is_strict_superset = partial(is_strict_superset, self.arr)

        update_wrapper(self.intersection, intersection)
        update_wrapper(self.union, union)
        update_wrapper(self.difference, difference)
        update_wrapper(self.symmetric_difference, symmetric_difference)
        update_wrapper(self.isdisjoint, isdisjoint)
        update_wrapper(self.issubset, issubset)
        update_wrapper(self.issuperset, issuperset)
        update_wrapper(self.isequal, isequal)
        update_wrapper(self.isunique, isunique)
        update_wrapper(self.is_strict_subset, is_strict_subset)
        update_wrapper(self.is_strict_superset, is_strict_superset)

    def __repr__(self):
        return f"<ArraySetOps for array {self.arr}>"

    def __str__(self):
        return self.__repr__()


def calculate_angle(point_a, point_b, point_c, degrees=True):
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

    v1_mag = np.linalg.norm(v1, axis=-1)
    v1_norm = v1 / v1_mag[:, :, np.newaxis]

    v2_mag = np.linalg.norm(v2, axis=-1)
    v2_norm = v2 / v2_mag[:, :, np.newaxis]

    res = np.sum(v1_norm * v2_norm, axis=-1)

    if degrees:
        return np.degrees(np.arccos(res))

    return np.arccos(res)


def array_distance(arr1, arr2):
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
    distance_array = np.linalg.norm(arr1[:, np.newaxis, :] - arr2, axis=-1)

    return distance_array


def matching_indices(arr1, arr2):
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
    # TODO: simplify this (just return the indices)
    return np.where(idx)[0]


def non_matching_indices(arr1, arr2):
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

    # TODO: simplify this (just return the indices)
    return np.where(idx)[0]
    # simplify return above
    # return idx


# These are all set operations we want to implement in our classes
def asvoid(arr):
    """
    Reference: https://stackoverflow.com/questions/16216078/test-for-membership-in-a-2d-numpy-array
    Based on http://stackoverflow.com/a/16973510/190597 (Jaime, 2013-06)
    View the array as dtype np.void (bytes). The items along the last axis are
    viewed as one value. This allows comparisons to be performed on the entire row.
    """
    arr = np.ascontiguousarray(arr)
    if np.issubdtype(arr.dtype, np.floating):
        """Care needs to be taken here since
        np.array([-0.]).view(np.void) != np.array([0.]).view(np.void)
        Adding 0. converts -0. to 0.
        """
        arr += 0.0
    return arr.view(np.dtype((np.void, arr.dtype.itemsize * arr.shape[-1])))


def intersection(arr1, arr2, assume_unique=False):
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

    arr1 = asvoid(arr1)
    arr2 = asvoid(arr2)
    return np.in1d(arr1, arr2, assume_unique)


def difference(arr1, arr2, assume_unique=False):
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

    arr1 = asvoid(arr1)
    arr2 = asvoid(arr2)

    return np.in1d(arr1, arr2, assume_unique, invert=True)


def symmetric_difference(arr1, arr2, assume_unique=False):
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


def union(arr1, arr2):
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


def isdisjoint(arr1, arr2):
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


def issubset(arr1, arr2):
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

    return bool(np.all(intersection(arr1, arr2)))


def issuperset(arr1, arr2):
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


def isequal(arr1, arr2):
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


def isunique(arr):
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

    return arr.shape[0] == np.unique(arr, axis=0).shape[0]


def is_strict_subset(arr1, arr2):
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


def is_strict_superset(arr1, arr2):
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
