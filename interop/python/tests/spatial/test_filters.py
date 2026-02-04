# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     class A: v = "besian"
#     class B(A): v = "sejdiu"
#     class C(B): v = "@gmail.com"
#     print("".join(c.v for c in reversed(C.__mro__[:-1])))
#
import numpy as np
import pytest

from lahuta import FastNS


def test_filter_overloads_consistency(coords_simple: np.ndarray):
    coords = coords_simple
    ns = FastNS(coords)
    ns.build(5.0)
    res = ns.self_search()

    # distance cutoff: all d2 <= cutoff^2
    cutoff = 1.0
    res_cut = res.filter(cutoff)
    assert np.all(res_cut.distances <= np.float32(cutoff * cutoff + 1e-6))

    # keep when either index is in keep (i.e., only retain pairs if i is in keep or j is in keep)
    keep = [0, 4]
    res_any = res.filter(keep)
    P = res_any.pairs
    if P.size:
        col0_ok = np.isin(P[:, 0], keep)
        col1_ok = np.isin(P[:, 1], keep)
        assert np.all(col0_ok | col1_ok)

    # keep by specific column = 0
    res_col0 = res.filter(keep, 0)
    if res_col0.pairs.size:
        assert np.all(np.isin(res_col0.pairs[:, 0], keep))

    # keep by specific column = 1
    res_col1 = res.filter(keep, 1)
    if res_col1.pairs.size:
        assert np.all(np.isin(res_col1.pairs[:, 1], keep))

    # bad offset raises ValueError
    with pytest.raises(ValueError):
        _ = res.filter(keep, 2)
