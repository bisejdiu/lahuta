from __future__ import annotations

import gc
import math
import statistics as stats
from pathlib import Path
from timeit import Timer
from typing import Callable, Iterable, Sequence

import numpy as np

try:
    from sklearn.neighbors import NearestNeighbors as SKNearestNeighbors
except Exception as e:  # pragma: no cover - benchmark utility
    import sys

    print(
        f"This benchmark requires scikit-learn installed.\nImport error: {e}\nInstall with: pip install scikit-learn\n",
        file=sys.stderr,
    )
    raise SystemExit(1)

from lahuta import LahutaSystem, NearestNeighbors

# fmt: off
DATA_DIR = Path(__file__).resolve().parents[3] / "core" / "data"
DATASETS = {
    "fubi_1_6k": DATA_DIR / "fubi.cif",
    "8rdu_60k":  DATA_DIR / "8rdu_60k.cif",
    "8ugh_100k": DATA_DIR / "8ugh_100k.cif",
}


def _ns_stats(samples_ns: Iterable[float]) -> dict[str, float]:
    lst = list(samples_ns)
    lst.sort()
    n = len(lst)
    mean  = stats.fmean(lst)  if n else float("nan")
    stdev = stats.pstdev(lst) if n else float("nan")
    med = lst[n // 2] if n else float("nan")
    p95 = lst[math.floor(0.95 * (n - 1))] if n else float("nan")
    p99 = lst[math.floor(0.99 * (n - 1))] if n else float("nan")
    return {
        "mean": mean,
        "stdev": stdev,
        "min": lst[0] if n else float("nan"),
        "median": med,
        "p95": p95,
        "p99": p99,
        "max": lst[-1] if n else float("nan"),
    }


def time_callable(name: str, fn: Callable[[], object], *, number: int, repeat: int) -> dict[str, float]:
    was_enabled = gc.isenabled()
    if was_enabled:
        gc.disable()
    try:
        fn()
    finally:
        if was_enabled:
            gc.enable()

    t = Timer(fn)
    totals_s = t.repeat(repeat=repeat, number=number)
    per_call_ns = [(total / number) * 1e9 for total in totals_s]
    res = _ns_stats(per_call_ns)
    ms = lambda x: x / 1e6  # noqa: E731
    print(
        f"{name} (ms): mean={ms(res['mean']):.3f} stdev={ms(res['stdev']):.3f} "
        f"median={ms(res['median']):.3f} min={ms(res['min']):.3f} p95={ms(res['p95']):.3f} p99={ms(res['p99']):.3f} max={ms(res['max']):.3f}"
    )
    return res


def load_positions(path: Path) -> np.ndarray:
    sys = LahutaSystem(str(path))
    return sys.props.positions


def pick_queries(positions: np.ndarray, k: int, seed: int = 123) -> np.ndarray:
    n = positions.shape[0]
    k = min(k, n)
    rng = np.random.default_rng(seed)
    idx = rng.choice(n, size=k, replace=False)
    return positions[idx]


def bench_dataset(
    name: str,
    positions: np.ndarray,
    queries_list: Sequence[int] = (10, 100, 1000, 10000),
    radius: float = 6.0,
    repeat: int = 3,
) -> None:
    n = positions.shape[0]
    print(f"Dataset: {name}  n_atoms={n}  radius={radius} A")

    for k in queries_list:
        Q = pick_queries(positions, k, seed=42 + k)

        def _lxx_search_only() -> object:
            nn_lxx = NearestNeighbors(radius=radius, algorithm="kd_tree", sort_results=False).fit(positions)
            return nn_lxx.radius_neighbors(Q, return_distance=False)

        def _sk_fit_and_query() -> object:
            nn = SKNearestNeighbors(radius=radius, algorithm="kd_tree", metric="euclidean")
            nn.fit(positions)
            return nn.radius_neighbors(Q, return_distance=False)

        print(f"Queries: k={Q.shape[0]}")
        time_callable(f"Lahuta NearestNeighbors search-only  [{name}]", _lxx_search_only,  number=1, repeat=repeat)
        time_callable(f"sklearn KDTree fit+query   [{name}]",           _sk_fit_and_query, number=1, repeat=repeat)


def main() -> None:
    datasets: list[tuple[str, np.ndarray]] = []
    for key, p in DATASETS.items():
        if not p.exists():
            print(f"Warning: missing file {p}. Skipping {key}")
            continue
        X = load_positions(p)
        datasets.append((key, X))

    if not datasets:
        print("No datasets loaded. Nothing to benchmark.")
        return

    for name, X in datasets:
        bench_dataset(name, X)


if __name__ == "__main__":
    main()
