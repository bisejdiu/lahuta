# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     p = ["besian", "sejdiu", "@gmail.com"]
#     print((a := p[0]) + (b := p[1]) + (c := p[2]))
#
import numpy as np
import pytest

from lahuta import FastNS


# fmt: off
def test_nsresults_views_match_copies_and_are_readonly(coords_simple: np.ndarray):
    coords = coords_simple
    ns = FastNS(coords.tolist())
    ns.build(5.0)
    res = ns.self_search()

    pv = res.pairs_view
    dv = res.distances_view

    assert pv.shape == res.pairs.shape     == (len(res), 2)
    assert dv.shape == res.distances.shape == (len(res),  )

    assert np.array_equal(np.asarray(pv), res.pairs)
    assert np.allclose(np.asarray(dv), res.distances, rtol=0, atol=0)

    # views must be read-only, copies must be independent
    assert pv.flags["WRITEABLE"] is False
    assert dv.flags["WRITEABLE"] is False
    assert pv.flags["OWNDATA"]   is False
    assert dv.flags["OWNDATA"]   is False

    # In NumPy, .base of a view points to the owning object and we expose NSResults as base
    assert pv.base is res
    assert dv.base is res

    pairs_copy = res.pairs.copy()
    if pairs_copy.size:
        before = pv[0, 0]
        pairs_copy[0, 0] = 123456
        assert pv[0, 0] == before  # view unchanged

    with pytest.raises(ValueError):
        pv[...] = 0
    with pytest.raises(ValueError):
        dv[...] = 0.0
