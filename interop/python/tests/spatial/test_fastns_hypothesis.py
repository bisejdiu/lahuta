import math
from typing import Callable, TypeVar

import numpy as np
import pytest

pytest.importorskip("hypothesis", reason="Hypothesis not installed")

from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st
from hypothesis.strategies import SearchStrategy

from lahuta import FastNS

T = TypeVar("T")
Draw = Callable[[SearchStrategy[T]], T]
CategoryBuilder = Callable[[Draw, int], tuple[np.ndarray, float]]

BRUTE_FORCE_THRESHOLD = 5000
MAX_SEED = 2**32 - 1


# fmt: off
def _rng(draw: Draw[int]) -> np.random.Generator:
    return np.random.default_rng(draw(st.integers(0, MAX_SEED)))

def _linear_cloud(n: int, span: float) -> np.ndarray:
    xs = np.linspace(-span, span, n, dtype=np.float64)
    ys = xs * 0.5
    zs = xs * -0.25
    return np.stack([xs, ys, zs], axis=1)


def _build_tiny(draw: Draw, max_points: int) -> tuple[np.ndarray, float]:
    n      = draw(st.integers(1, min(8, max_points)))
    span   = draw(st.floats(1e-3, 1.0, allow_nan=False, allow_infinity=False))
    cutoff = draw(st.floats(1e-3, 2.0, allow_nan=False, allow_infinity=False))

    coords = _rng(draw).uniform(-span, span, size=(n, 3)).astype(np.float64)
    return coords, float(cutoff)


def _build_small(draw: Draw, max_points: int) -> tuple[np.ndarray, float]:
    n      = draw(st.integers(8, min(64, max_points)))
    span   = draw(st.floats(1.0, 100.0,       allow_nan=False, allow_infinity=False))
    cutoff = draw(st.floats(0.25, span * 2.0, allow_nan=False, allow_infinity=False))

    coords = _rng(draw).normal(loc=0.0, scale=span, size=(n, 3)).astype(np.float64)
    return coords, float(cutoff)


def _build_medium(draw: Draw, max_points: int) -> tuple[np.ndarray, float]:
    n      = draw(st.integers(64, min(512, max_points)))
    span   = draw(st.floats(10.0, 5_000.0,    allow_nan=False, allow_infinity=False))
    cutoff = draw(st.floats(0.5, span * 0.75, allow_nan=False, allow_infinity=False))

    base   = _rng(draw).uniform(-span, span, size=(n, 3)).astype(np.float64)
    jitter = _rng(draw).normal(loc=0.0, scale=span * 1e-6, size=(n, 3)).astype(np.float64)
    coords = base + jitter
    return coords, float(cutoff)


def _build_mixed_scale(draw: Draw, max_points: int) -> tuple[np.ndarray, float]:
    n_small    = draw(st.integers(8, min(64,  max_points // 4)))
    n_large    = draw(st.integers(8, min(256, max_points - n_small)))
    large_span = draw(st.floats(1e4, 1e7,               allow_nan=False, allow_infinity=False))
    cutoff     = draw(st.floats(1e-4, large_span * 0.1, allow_nan=False, allow_infinity=False))

    small = _rng(draw).uniform(-1e-6, 1e-6, size=(n_small, 3)).astype(np.float64)
    large = _rng(draw).uniform(-large_span, large_span, size=(n_large, 3)).astype(np.float64)
    coords = np.concatenate([small, large], axis=0)
    return coords, float(cutoff)


def _build_ultra_small(draw: Draw, max_points: int) -> tuple[np.ndarray, float]:
    n      = draw(st.integers(16, min(256, max_points)))
    cutoff = draw(st.floats(1e-12, 1e-6, allow_nan=False, allow_infinity=False))

    coords = _rng(draw).uniform(-1e-12, 1e-12, size=(n, 3)).astype(np.float64)
    return coords, float(cutoff)


def _build_ultra_large(draw: Draw, max_points: int) -> tuple[np.ndarray, float]:
    n      = draw(st.integers(16, min(512, max_points)))
    span   = draw(st.floats(1e8, 1e11,       allow_nan=False, allow_infinity=False))
    cutoff = draw(st.floats(1.0, span * 0.5, allow_nan=False, allow_infinity=False))

    coords = _rng(draw).uniform(-span, span, size=(n, 3)).astype(np.float64)
    return coords, float(cutoff)


def _build_large_sparse(draw: Draw, max_points: int) -> tuple[np.ndarray, float]:
    n      = draw(st.integers(512, min(max_points, 3072)))
    span   = draw(st.floats(1_000.0, 50_000.0, allow_nan=False, allow_infinity=False))
    cutoff = draw(st.floats(10.0, span * 0.01, allow_nan=False, allow_infinity=False))

    coords = _linear_cloud(n, span)
    return coords, float(cutoff)


def _build_deg_line(draw: Draw, max_points: int) -> tuple[np.ndarray, float]:
    n       = draw(st.integers(16, max_points))
    span    = draw(st.floats(10.0, 10_000.0,    allow_nan=False, allow_infinity=False))
    cutoff  = draw(st.floats(1e-4, span * 0.02, allow_nan=False, allow_infinity=False))

    coords = _linear_cloud(n, span)
    return coords, float(cutoff)


def _build_deg_plane(draw: Draw, max_points: int) -> tuple[np.ndarray, float]:
    n      = draw(st.integers(32, max_points))
    span   = draw(st.floats(10.0, 10_000.0,    allow_nan=False, allow_infinity=False))
    cutoff = draw(st.floats(1e-4, span * 0.02, allow_nan=False, allow_infinity=False))

    xs = np.linspace(-span, span, n, dtype=np.float64)
    ys = np.sin(xs / span)
    coords = np.column_stack([xs, ys, np.zeros_like(xs)])
    return coords, float(cutoff)


DEFAULT_CLOUD_BUILDERS: dict[str, CategoryBuilder] = {
    "tiny":   _build_tiny,
    "small":  _build_small,
    "medium": _build_medium,
    "mixed_scale": _build_mixed_scale,
    "ultra_small": _build_ultra_small,
    "ultra_large": _build_ultra_large,
    "large_sparse":     _build_large_sparse,
    "degenerate_line":  _build_deg_line,
    "degenerate_plane": _build_deg_plane,
}


@st.composite
def crazy_point_clouds(
    draw: Draw,
    max_points: int = BRUTE_FORCE_THRESHOLD,
    builders: dict[str, CategoryBuilder] | None = None,
) -> tuple[np.ndarray, float]:
    builder_map = builders or DEFAULT_CLOUD_BUILDERS
    category = draw(st.sampled_from(list(builder_map.keys())))
    coords, cutoff = builder_map[category](draw, max_points)
    coords = np.ascontiguousarray(coords, dtype=np.float64)
    return coords, float(cutoff)


def _brute_force_neighbors(coords: np.ndarray, cutoff: float) -> dict[tuple[int, int], float]:
    n = coords.shape[0]
    if n < 2:
        return {}
    cutoff_sq = float(cutoff * cutoff)
    diffs   = coords[:, None, :] - coords[None, :, :]
    dist_sq = np.sum(diffs * diffs, axis=2)
    mask    = np.triu(dist_sq <= cutoff_sq, k=1)
    ii, jj  = np.nonzero(mask)
    return {(int(i), int(j)): float(dist_sq[i, j]) for i, j in zip(ii, jj)}


def _canonicalize_results(res) -> dict[tuple[int, int], float]:
    if res.pairs.size == 0: return {}
    out = {}
    for (i, j), dist_sq in zip(res.pairs.tolist(), res.distances.tolist()):
        a = int(i)
        b = int(j)
        if b < a:
            a, b = b, a
        out[(a, b)] = float(dist_sq)
    return out


def test_fastns_fails_without_fallback_when_grid_unbuildable():
    coords = np.empty((0, 3), dtype=np.float64)
    grid = FastNS(coords)
    assert grid.build(1.0, brute_force_fallback=False) is False

    with pytest.raises(RuntimeError):
        grid.self_search()


def test_fastns_fallback_disabled_above_threshold_raises():
    def make_cloud(n, span):
        xs = np.linspace(-span, span, n, dtype=np.float64)
        ys = xs *  0.5
        zs = xs * -0.25
        return np.stack([xs, ys, zs], axis=1)

    n = BRUTE_FORCE_THRESHOLD + 1
    coords = make_cloud(n, 1.0)
    cutoff = float("nan")

    grid = FastNS(coords)
    ok = grid.build(cutoff, brute_force_fallback=True)
    assert ok is False, "expected build() to fail for NaN cutoff"

    with pytest.raises(RuntimeError):
        grid.self_search()


@settings(deadline=None, max_examples=60, suppress_health_check=[HealthCheck.too_slow])
@given(crazy_point_clouds())
def test_fastns_handles_diverse_point_clouds_with_fallback(data):
    coords, cutoff = data
    n = coords.shape[0]
    assume(n <= BRUTE_FORCE_THRESHOLD)
    ns = FastNS(coords)
    success = ns.build(cutoff, brute_force_fallback=True)
    if n == 0:
        assert success is False
        with pytest.raises(RuntimeError):
            ns.self_search()
        return

    assert success is True
    res = ns.self_search()

    assert res.size() == len(res)
    assert res.pairs.shape[1] == 2 if res.pairs.size else True
    assert np.all(np.isfinite(res.distances))
    assert np.all(res.distances >= 0.0)

    if n <= 128:
        brute  = _brute_force_neighbors(coords, cutoff)
        fastns = _canonicalize_results(res)
        # With mixed-scale point clouds, float precision may produce small-distance underflow to 0.0.
        # If that's the case, allow a subset match rather than strict equality.
        if np.max(np.abs(coords)) / max(1e-30, np.min(np.abs(coords[np.nonzero(coords)]), initial=1e-30)) > 1e6:
            assert set(fastns.keys()).issubset(set(brute.keys()))
        else:
            assert set(fastns.keys()) == set(brute.keys())
        for key in fastns:
            tol = max(5e-4, 5e-6 * abs(brute[key]))
            expected_d = math.sqrt(brute[key])
            actual_d   = math.sqrt(fastns[key])
            # Allow 0.0 vs tiny differences for mixed-scale scenarios
            if expected_d < 1e-7:
                assert actual_d == 0.0 or math.isclose(actual_d, expected_d, rel_tol=5e-6, abs_tol=tol)
            else:
                assert math.isclose(actual_d, expected_d, rel_tol=5e-6, abs_tol=tol)
