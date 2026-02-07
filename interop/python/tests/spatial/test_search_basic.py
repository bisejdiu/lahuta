# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     f = lambda g: lambda a: lambda b: lambda c: g(a, b, c)
#     g = lambda a, b, c: a + b + c
#     print(f(g)("besian")("sejdiu")("@gmail.com"))
#
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


def test_kdindex_build_view_accepts_float64(coords_simple: np.ndarray) -> None:
    coords = coords_simple
    assert coords.dtype == np.float64
    assert coords.flags.c_contiguous

    kd = KDIndex()
    assert kd.build_view(coords)
    assert kd.ready


def test_kdindex_radius_search_grouped_matches_nsresults(coords_simple: np.ndarray) -> None:
    coords = coords_simple
    kd = KDIndex()
    assert kd.build(coords)

    def expected_grouped(
        pairs: np.ndarray,
        d2: np.ndarray,
        *,
        return_distance: bool,
        sort_results: bool,
    ) -> list[np.ndarray] | tuple[list[np.ndarray], list[np.ndarray]]:
        n_queries = coords.shape[0]
        idx_lists: list[list[int]] = [[] for _ in range(n_queries)]
        dist_lists: list[list[float]] = [[] for _ in range(n_queries)]

        for (i, j), v in zip(pairs.tolist(), d2.tolist()):
            idx_lists[i].append(int(j))
            dist_lists[i].append(math.sqrt(max(0.0, float(v))))

        indices: list[np.ndarray] = []
        distances: list[np.ndarray] = []
        for idxs, dsts in zip(idx_lists, dist_lists):
            if sort_results:
                if return_distance:
                    order = np.argsort(np.asarray(dsts, dtype=np.float64))
                    idxs = [idxs[int(o)] for o in order]
                    dsts = [dsts[int(o)] for o in order]
                else:
                    idxs = sorted(idxs)
            indices.append(np.asarray(idxs, dtype=np.int32))
            distances.append(np.asarray(dsts, dtype=np.float64))

        return (distances, indices) if return_distance else indices

    flat = kd.radius_search(coords, 1.1)
    expected = expected_grouped(
        flat.pairs,
        flat.distances,
        return_distance=False,
        sort_results=False,
    )
    grouped = kd.radius_search(coords, 1.1, grouped=True, sort_results=False)
    assert isinstance(grouped, list)
    assert len(grouped) == len(expected) == coords.shape[0]
    for g, e in zip(grouped, expected):
        np.testing.assert_array_equal(g, e)

    expected_distances, expected_indices = expected_grouped(
        flat.pairs,
        flat.distances,
        return_distance=True,
        sort_results=False,
    )
    grouped_distances, grouped_indices = kd.radius_search(
        coords, 1.1, grouped=True, return_distance=True, sort_results=False
    )
    assert len(grouped_distances) == len(grouped_indices) == coords.shape[0]
    assert len(expected_distances) == len(expected_indices) == coords.shape[0]
    for gd, gi, ed, ei in zip(grouped_distances, grouped_indices, expected_distances, expected_indices):
        np.testing.assert_array_equal(gi, ei)
        np.testing.assert_allclose(gd, ed, rtol=1e-6, atol=1e-6)

    sorted_distances, sorted_indices = kd.radius_search(
        coords, 1.1, grouped=True, return_distance=True, sort_results=True
    )
    for sd, si in zip(sorted_distances, sorted_indices):
        assert sd.shape == si.shape
        if sd.size > 1:
            assert np.all(sd[:-1] <= sd[1:])


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
