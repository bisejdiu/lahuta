import math

import numpy as np
import pytest

pytest.importorskip("hypothesis", reason="Hypothesis not installed")

from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra import numpy as hnp

from lahuta import FastNS, NSResults


# fmt: off
def _pairs_to_dict(pairs: np.ndarray, d2: np.ndarray) -> dict[tuple[int, int], float]:
    m: dict[tuple[int, int], float] = {}
    for (i, j), dd in zip(pairs.tolist(), d2.tolist()):
        m[(int(i)), int(j)] = float(dd)
    return m


def _assert_equivalent(a: NSResults, b: NSResults, rtol=1e-6, atol=1e-6):
    amap = _pairs_to_dict(a.pairs, a.distances)
    bmap = _pairs_to_dict(b.pairs, b.distances)
    assert set(amap.keys()) == set(bmap.keys())
    for k in amap:
        assert math.isclose(amap[k], bmap[k], rel_tol=rtol, abs_tol=atol)


# Generate random coordinate arrays of shape (n,3), n in [0, 50], values bounded
coords_arrays = hnp.arrays(
    dtype=np.float64,
    shape=hnp.array_shapes(min_dims=2, max_dims=2, min_side=0, max_side=50).map(lambda s: (s[0], 3)),
    elements=st.floats(allow_nan=False, allow_infinity=False, width=64, min_value=-5.0, max_value=5.0),
)


@settings(deadline=None, max_examples=50)
@given(coords=coords_arrays, cutoff=st.floats(min_value=0.5, max_value=3.0))
def test_constructors_equivalence_and_view_copy_consistency(coords: np.ndarray, cutoff: float):
    from lahuta import rdkit as lrdkit

    ns_np   = FastNS(coords)
    ns_list = FastNS(coords.tolist())
    ns_pts  = FastNS([lrdkit.Point3D(*row) for row in coords.tolist()])

    ok_np   = ns_np  .build(cutoff)
    ok_list = ns_list.build(cutoff)
    ok_pts  = ns_pts .build(cutoff)

    assert ok_np == ok_list == ok_pts
    if not ok_np: return

    res_np   = ns_np  .self_search()
    res_list = ns_list.self_search()
    res_pts  = ns_pts .self_search()

    _assert_equivalent(res_np, res_list)
    _assert_equivalent(res_np, res_pts)

    res_round = NSResults(res_np.pairs, res_np.distances)
    _assert_equivalent(res_np, res_round)

    pv = res_np.pairs_view
    dv = res_np.distances_view
    assert pv.shape == res_np.pairs.shape
    assert dv.shape == res_np.distances.shape
    assert np.array_equal(np.asarray(pv), res_np.pairs)
    assert np.allclose(np.asarray(dv), res_np.distances, rtol=0, atol=0)
    assert pv.flags["WRITEABLE"] is False
    assert dv.flags["WRITEABLE"] is False

    if len(res_np) > 0:
        sqrt_calc = np.sqrt(np.maximum(0.0, res_np.distances.astype(np.float64))).astype(np.float32)
        assert np.allclose(res_np.get_sqrt_distances(), sqrt_calc, rtol=1e-6, atol=1e-6)

    # Filter invariants
    out_cut = res_np.filter(cutoff)
    assert np.all(out_cut.distances <= np.float32(cutoff * cutoff + 1e-6))

    n = coords.shape[0]
    if n > 0 and len(res_np) > 0:
        S = set(np.random.default_rng(0).choice(n, size=max(1, n // 3), replace=False).tolist())
        out_any = res_np.filter(sorted(S))
        if len(out_any) > 0:
            P = out_any.pairs
            assert np.all(np.isin(P[:, 0], list(S)) | np.isin(P[:, 1], list(S)))

        out_col0 = res_np.filter(sorted(S), 0)
        if len(out_col0) > 0:
            assert np.all(np.isin(out_col0.pairs[:, 0], list(S)))

        out_col1 = res_np.filter(sorted(S), 1)
        if len(out_col1) > 0:
            assert np.all(np.isin(out_col1.pairs[:, 1], list(S)))


@settings(deadline=None, max_examples=25)
@given(
    pairs = hnp.arrays(dtype=np.int32,   shape=hnp.array_shapes(min_dims=2, max_dims=2, min_side=0, max_side=40).map(lambda s: (s[0], 2))),
    d2    = hnp.arrays(dtype=np.float32, shape=hnp.array_shapes(min_dims=1, max_dims=1, min_side=0, max_side=40).map(lambda s: (s[0],  ))),
)
def test_nsresults_numpy_ctor_shape_checks_and_sqrt_clamp(pairs: np.ndarray, d2: np.ndarray):
    if pairs.shape[0] != d2.shape[0]:
        with pytest.raises(ValueError):
            NSResults(pairs, d2)
        return

    res = NSResults(pairs, d2)
    assert res.pairs.shape     == pairs.shape
    assert res.distances.shape == d2.shape

    if res.distances.size > 0:
        d2_mod = res.distances.copy()
        d2_mod[0] = np.float32(-1e-10)  # injectng a tiny negative value to test clamping
        res2 = NSResults(res.pairs, d2_mod)
        sqrt2 = res2.get_sqrt_distances()
        assert sqrt2[0] == np.float32(0.0)
