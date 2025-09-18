import math

import numpy as np
import pytest

from lahuta import FastNS, KDIndex, NSResults


# fmt: off
def test_fastns_build_and_self_search_basic(coords_simple: np.ndarray):
    coords = coords_simple
    ns = FastNS(coords)
    ok = ns.build(1.1)
    assert ok is True
    res = ns.self_search()
    assert isinstance(res, NSResults)

    # sanity checks
    assert len(res) == res.size()
    assert res.pairs.shape[1] == 2
    assert res.distances.ndim == 1

    # types
    assert np.issubdtype(res.pairs.dtype, np.integer)
    assert res.distances.dtype == np.float32


def test_iteration_protocol_matches_pairs_and_distances(coords_simple: np.ndarray):
    coords = coords_simple
    ns = FastNS(coords)
    ns.build(5.0)
    res = ns.self_search()

    seen = []
    for (i, j), d2 in res:
        seen.append(((i, j), float(d2)))
    assert len(seen) == len(res)

    P  = res.pairs
    D2 = res.distances
    sqrt_all = res.get_sqrt_distances()
    for idx, ((i, j), d2) in enumerate(seen):
        assert i == int(P[idx, 0]) and j == int(P[idx, 1])
        assert math.isclose(d2, float(D2[idx]), rel_tol=0, abs_tol=0)
        assert math.isclose(math.sqrt(max(0.0, d2)), float(sqrt_all[idx]), rel_tol=1e-6, abs_tol=1e-6)


def test_search_numpy_overload_equivalence():
    rng = np.random.default_rng(0)
    coords = rng.normal(size=(50, 3)).astype(np.float64)

    ns = FastNS(coords)
    assert ns.build(2.0)

    res_self = ns.self_search()
    kd = KDIndex()
    assert kd.build(coords.tolist())
    res_np   = kd.radius_search(coords, 2.0)

    # self_search returns unique pairs among stored coords, search(coords) returns
    # (query_idx, stored_idx), including identity pairs (i, i) and mirrored pairs.
    # Canonicalize: drop identity and map both to unordered pairs.
    def canon_pairs(pairs: np.ndarray, d2: np.ndarray) -> dict[tuple[int, int], float]:
        out: dict[tuple[int, int], float] = {}
        for (i, j), v in zip(pairs.tolist(), d2.tolist()):
            if int(i) == int(j): continue

            a, b = (int(i), int(j))
            key  = (a, b) if a < b else (b, a)

            # distances should be consistent for mirrored pairs, last write wins
            out[key] = float(v)
        return out

    a = canon_pairs(res_self.pairs, res_self.distances)
    b = canon_pairs(res_np.pairs,   res_np.distances)
    assert set(a.keys()) == set(b.keys())
    for k in a:
        assert math.isclose(a[k], b[k], rel_tol=1e-6, abs_tol=1e-6)


def test_search_numpy_bad_shape_errors():
    rng = np.random.default_rng(1)
    coords = rng.normal(size=(10, 3)).astype(np.float64)
    ns = FastNS(coords)
    assert ns.build(1.0)

    bad = np.ones((10, 2), dtype=np.float64)

    kd = KDIndex()
    assert kd.build(coords)
    with pytest.raises((ValueError, TypeError)):
        kd.radius_search(bad, 1.0)
