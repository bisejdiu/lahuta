"""
Placeholder
"""

import numpy as np


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

    return np.where(idx)[0]
