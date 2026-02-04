# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     @functools.wraps(lambda: None, updated=())
#     def f():
#         return "besian" + "sejdiu" + "@gmail.com"
#     print(f())
#
from typing import Any

import numpy as np
import pytest

from lahuta import metrics, NearestNeighbors

_spdist: Any
_sk_pairwise_distances: Any
_SKNearestNeighbors: Any

try:
    import scipy.spatial.distance as _spdist

    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

try:
    from sklearn.metrics import pairwise_distances as _sk_pairwise_distances
    from sklearn.neighbors import NearestNeighbors as _SKNearestNeighbors

    HAVE_SKLEARN = True
except Exception:
    HAVE_SKLEARN = False


# At some point we may need to tighten up how we handle these tolerances.
# See: core/tests/test_utils/fp_test_utils.hpp
TOL_R = 5e-7
TOL_A = 1e-12


def _cdist_numpy(X: np.ndarray, Y: np.ndarray, squared: bool = False) -> np.ndarray:
    """Numpy implementation for cross-distance computation."""
    diff = X[:, None, :] - Y[None, :, :]
    d2 = np.sum(diff * diff, axis=-1)
    return d2 if squared else np.sqrt(d2)


def _pdist_numpy(X: np.ndarray, squared: bool = False) -> np.ndarray:
    """Numpy implementation for pairwise distance computation."""
    D = _cdist_numpy(X, X, squared=squared)
    n = D.shape[0]
    iu = np.triu_indices(n, k=1)
    return D[iu]


def _pairwise_numpy(X: np.ndarray, Y: np.ndarray | None, squared: bool) -> np.ndarray:
    """Numpy implementation for pairwise distance matrix."""
    return _cdist_numpy(X, X if Y is None else Y, squared=squared)


def ref_cdist(X: np.ndarray, Y: np.ndarray, squared: bool) -> np.ndarray:
    """Reference implementation for cross-distance computation."""
    if HAVE_SCIPY:
        D = _spdist.cdist(X, Y, metric="euclidean")
        return D * D if squared else D
    return _cdist_numpy(X, Y, squared=squared)


def ref_pdist(X: np.ndarray, squared: bool) -> np.ndarray:
    """Reference implementation for pairwise distance computation."""
    if HAVE_SCIPY:
        d = _spdist.pdist(X, metric="euclidean")
        return d * d if squared else d
    return _pdist_numpy(X, squared=squared)


def ref_pairwise(X: np.ndarray, Y: np.ndarray | None, squared: bool) -> np.ndarray:
    """Reference implementation for pairwise distance matrix."""
    if HAVE_SKLEARN:
        D = _sk_pairwise_distances(X, X if Y is None else Y, metric="euclidean")
        return D * D if squared else D
    return _pairwise_numpy(X, Y, squared=squared)


def ref_radius_neighbors(X: np.ndarray, radius: float, Y: np.ndarray | None):
    """Reference implementation for radius neighbors search."""
    # returns (dists_list, idxs_list)
    if HAVE_SKLEARN:
        A = X if Y is None else Y
        nn = _SKNearestNeighbors(radius=radius, algorithm="brute", metric="euclidean")
        nn.fit(A)
        dists, idxs = nn.radius_neighbors(X, return_distance=True)
        out_d, out_i = [], []
        for i in range(X.shape[0]):
            di = dists[i]
            ii = idxs[i]
            if Y is None:
                mask = ii != i
                di = di[mask]
                ii = ii[mask]
            if di.size > 1:
                order = np.argsort(di)
                di = di[order]
                ii = ii[order]
            out_d.append(di)
            out_i.append(ii)
        return out_d, out_i
    # Numpy fallback
    B = X if Y is None else Y
    out_i, out_d = [], []
    for i in range(X.shape[0]):
        diff = B - X[i]
        d = np.sqrt(np.sum(diff * diff, axis=1))  # slightly faster than np.linalg.norm
        if Y is None:
            mask = (d <= radius) & (np.arange(B.shape[0]) != i)
        else:
            mask = d <= radius
        sel = np.where(mask)[0]
        di = d[sel]
        order = np.argsort(di)
        out_i.append(sel[order].astype(np.int32))
        out_d.append(di[order].astype(np.float64))
    return out_d, out_i


@pytest.mark.parametrize("squared", [False, True])
def test_cdist_matches_reference(squared: bool) -> None:
    """Test that cdist matches reference implementation."""
    rng = np.random.default_rng(42)
    X = rng.standard_normal((13, 3)).astype(np.float64)
    Y = rng.standard_normal((7, 3)).astype(np.float64)

    D_ref = ref_cdist(X, Y, squared=squared)
    D = metrics.cdist(X, Y, squared=squared)
    assert D.shape == D_ref.shape
    assert np.allclose(D, D_ref, rtol=TOL_R, atol=TOL_A)


def test_pdist_matches_reference() -> None:
    """Test that pdist matches reference implementation."""
    rng = np.random.default_rng(0)
    X = rng.random((10, 3), dtype=np.float64)
    for squared in (False, True):
        d_ref = ref_pdist(X, squared=squared)
        d = metrics.pdist(X, squared=squared)
        assert d.shape == d_ref.shape
        assert np.allclose(d, d_ref, rtol=TOL_R, atol=TOL_A)


def test_pairwise_distances_matches_reference() -> None:
    """Test that pairwise_distances matches reference implementation."""
    rng = np.random.default_rng(123)
    X = rng.standard_normal((9, 3)).astype(np.float64)
    Y = rng.standard_normal((5, 3)).astype(np.float64)
    for squared in (False, True):
        D_ref = ref_pairwise(X, None, squared=squared)
        D = metrics.pairwise_distances(X, squared=squared)
        assert np.allclose(D, D_ref, rtol=TOL_R, atol=TOL_A)

        C_ref = ref_pairwise(X, Y, squared=squared)
        C = metrics.pairwise_distances(X, Y, squared=squared)
        assert np.allclose(C, C_ref, rtol=TOL_R, atol=TOL_A)


def test_radius_neighbors_matches_reference_self_and_cross() -> None:
    """Test that radius_neighbors matches reference implementation for self and cross queries."""
    rng = np.random.default_rng(7)
    X = rng.standard_normal((20, 3)).astype(np.float64)
    Y = rng.standard_normal((15, 3)).astype(np.float64)
    radius = 1.25

    # self via sklearn-like wrapper
    d_ref, i_ref = ref_radius_neighbors(X, radius, Y=None)
    nn_self = NearestNeighbors(radius=radius, algorithm="kd_tree", sort_results=True).fit(X)
    d_lst, i_lst = nn_self.radius_neighbors(return_distance=True)
    assert len(d_lst) == len(i_lst) == X.shape[0]
    for dr, ir, dl, il in zip(d_ref, i_ref, d_lst, i_lst):
        # Both are sorted ascending
        # Compare content
        np.testing.assert_array_equal(il, np.asarray(ir, dtype=np.int32))
        np.testing.assert_allclose(dl, np.asarray(dr, dtype=np.float64), rtol=TOL_R, atol=TOL_A)

    # cross via wrapper (KD index on Y)
    d_ref_c, i_ref_c = ref_radius_neighbors(X, radius, Y=Y)
    nn_cross = NearestNeighbors(radius=radius, algorithm="kd_tree", sort_results=True).fit(Y)
    d_lst_c, i_lst_c = nn_cross.radius_neighbors(X, return_distance=True)
    assert len(d_lst_c) == len(i_lst_c) == X.shape[0]
    for dr, ir, dl, il in zip(d_ref_c, i_ref_c, d_lst_c, i_lst_c):
        # Reference may not be sorted
        # Sort pairs by distance
        if len(dr) > 1:
            order = np.argsort(dr)
            dr = np.asarray(dr)[order]
            ir = np.asarray(ir)[order]
        np.testing.assert_array_equal(il, np.asarray(ir, dtype=np.int32))
        np.testing.assert_allclose(dl, np.asarray(dr, dtype=np.float64), rtol=TOL_R, atol=TOL_A)
