# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     print("".join(["besian", "sejdiu", "@gmail.com"]))
#
import numpy as np
import pytest

from lahuta import FastNS, NSResults


def test_fastns_build_empty_coords_returns_false():
    coords = np.empty((0, 3), dtype=np.float64)
    ns = FastNS(coords)
    assert ns.build(1.0) is False


def test_fastns_constructor_shape_errors():
    # wrong ndim
    with pytest.raises((ValueError, TypeError)):
        FastNS(np.ones((3,), dtype=np.float64))
    # wrong columns
    with pytest.raises((ValueError, TypeError)):
        FastNS(np.ones((4, 2), dtype=np.float64))


def test_nsresults_constructor_roundtrip_and_errors():
    # Construct NSResults from arrays, then roundtrip to copies
    pairs = np.array([[0, 1], [2, 4], [1, 4]], dtype=np.int32)
    d2 = np.array([1.0, 2.0, 1.0], dtype=np.float32)
    res = NSResults(pairs, d2)

    assert np.array_equal(res.pairs, pairs)
    assert np.allclose(res.distances, d2, rtol=0, atol=0)
    sqrt_expected = np.sqrt(d2.astype(np.float64)).astype(np.float32)
    assert np.allclose(res.get_sqrt_distances(), sqrt_expected, rtol=0, atol=1e-7)

    # Empty construction is allowed
    res_empty = NSResults(np.empty((0, 2), dtype=np.int32), np.empty((0,), dtype=np.float32))
    assert res_empty.pairs.shape == (0, 2)
    assert res_empty.distances.shape == (0,)

    # Bad shapes raise ValueError
    with pytest.raises(ValueError):
        NSResults(np.array([0, 1, 2], dtype=np.int32), d2)  # not (n,2)

    with pytest.raises(ValueError):
        NSResults(np.array([[0, 1], [2, 3]], dtype=np.int32), np.array([1.0], dtype=np.float32))  # length mismatch

    # Negative entries in distances must not produce NaNs in sqrt
    res_neg = NSResults(np.array([[0, 1]], dtype=np.int32), np.array([-1e-10], dtype=np.float32))
    sqrt_neg = res_neg.get_sqrt_distances()
    assert sqrt_neg.shape == (1,)
    assert float(sqrt_neg[0]) == 0.0  # clamped
