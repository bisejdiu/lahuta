# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     print(string.Template("$a$b@$c").substitute(a="besian", b="sejdiu", c="gmail.com"))
#
from __future__ import annotations

import gc
import math
import statistics as stats
from pathlib import Path
from timeit import Timer
from typing import Callable, Iterable

import numpy as np

try:
    from sklearn.neighbors import NearestNeighbors as SKNearestNeighbors
except Exception as e:
    import sys

    print(
        f"This benchmark requires scikit-learn installed.\nImport error: {e}\nInstall with: pip install scikit-learn\n",
        file=sys.stderr,
    )
    raise SystemExit(1)

from lahuta import FastNS, LahutaSystem, NearestNeighbors

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
    # Warmup outside measurement. GC disabled around timeit runs
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
    return sys.props.positions  # (n, 3) float64 copy


def bench_self(name: str, positions: np.ndarray, radius: float = 6.0, repeat: int = 1) -> None:
    n = positions.shape[0]
    print(f"Self-search: {name}  n_atoms={n}  radius={radius} A")

    def _lahuta_nn_self() -> object:
        nn = NearestNeighbors(radius=radius, algorithm="kd_tree", sort_results=False).fit(positions)
        return nn.radius_neighbors(return_distance=False)

    def _lahuta_ns_self() -> object:
        ns = FastNS(positions)
        ns.build(radius)
        return ns.self_search()

    def _sk_self() -> object:
        nn = SKNearestNeighbors(radius=radius, algorithm="kd_tree", metric="euclidean")
        nn.fit(positions)
        return nn.radius_neighbors(positions, return_distance=False)

    time_callable(f"Lahuta NearestNeighbors fit+radius_neighbors self [{name}]", _lahuta_nn_self, number=1, repeat=repeat)
    time_callable(f"Lahuta FastNS build+self_search [{name}]",                   _lahuta_ns_self, number=1, repeat=repeat)
    time_callable(f"sklearn KDTree fit+radius_neighbors self [{name}]",         _sk_self, number=1, repeat=repeat)

    lahuta_idxs = NearestNeighbors(radius=radius, algorithm="kd_tree").fit(positions).radius_neighbors(return_distance=False)
    lahuta_total_links = int(sum(map(len, lahuta_idxs)))
    lahuta_unique_pairs = lahuta_total_links // 2

    skl_nn = SKNearestNeighbors(radius=radius, algorithm="kd_tree", metric="euclidean").fit(positions)
    skl_idxs = skl_nn.radius_neighbors(positions, return_distance=False)
    total_links = int(sum(map(len, skl_idxs)))
    sk_unique_pairs = (total_links - n) // 2

    print(f"Pairs (unique, undirected): Lahuta={lahuta_unique_pairs:,}  sklearn~={sk_unique_pairs:,}")


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

    # repeat=1 because these are slow
    for name, X in datasets:
        bench_self(name, X, radius=6.0, repeat=1)


if __name__ == "__main__":
    main()
