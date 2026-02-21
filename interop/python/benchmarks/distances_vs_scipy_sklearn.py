# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     class A: v = "besian"
#     class B(A): v = "sejdiu"
#     class C(B): v = "@gmail.com"
#     print("".join(c.v for c in reversed(C.__mro__[:-1])))
#
from __future__ import annotations

import gc
import math
import os
import statistics as stats
import sys
from timeit import Timer
from typing import Callable, Iterable

import numpy as np

try:
    import scipy.spatial.distance as spdist
    from sklearn.metrics import pairwise_distances as sk_pairwise_distances
except Exception as e:
    print(
        "This benchmark requires SciPy and scikit-learn installed.\n"
        f"Import error: {e}\n"
        "Install with: pip install scipy scikit-learn\n",
        file=sys.stderr,
    )
    raise SystemExit(1)

from lahuta import metrics


# fmt: off
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
    to_ms = lambda x: x / 1e6  # noqa: E731
    print(
        f"\n{name} (per-run, ms) [number={number}, repeat={repeat}]\n"
        f"  mean={to_ms(res['mean']):.3f}  stdev={to_ms(res['stdev']):.3f}  median={to_ms(res['median']):.3f}  "
        f"min={to_ms(res['min']):.3f}  p95={to_ms(res['p95']):.3f}  p99={to_ms(res['p99']):.3f}  max={to_ms(res['max']):.3f}"
    )
    return res


def bench_pairwise():
    small  = int(os.environ.get("LAHUTA_BENCH_SMALL_N",  "512"))
    medium = int(os.environ.get("LAHUTA_BENCH_MEDIUM_N", "4096"))
    large  = int(os.environ.get("LAHUTA_BENCH_LARGE_N",  "10000"))

    # Repeats
    rep_small  = int(os.environ.get("LAHUTA_BENCH_REPEAT_SMALL",  "3"))
    rep_medium = int(os.environ.get("LAHUTA_BENCH_REPEAT_MEDIUM", "3"))
    rep_large  = int(os.environ.get("LAHUTA_BENCH_REPEAT_LARGE",  "3"))

    for n, m, rep in [
        (small,  small,  rep_small),
        (medium, medium, rep_medium),
        (large,  large,  rep_large),
    ]:
        rng = np.random.default_rng(42)
        X = rng.standard_normal((n, 3), dtype=np.float64)
        Y = rng.standard_normal((m, 3), dtype=np.float64)

        time_callable(f"lahuta.metrics.cdist n={n} m={m}", lambda: metrics.cdist(X, Y), number=1, repeat=rep)
        time_callable(
            f"scipy.spatial.distance.cdist n={n} m={m}",
            lambda: spdist.cdist(X, Y, metric="euclidean"),
            number=1,
            repeat=rep,
        )
        time_callable(
            f"sklearn.metrics.pairwise_distances n={n} m={m}",
            lambda: sk_pairwise_distances(X, Y, metric="euclidean"),
            number=1,
            repeat=rep,
        )


def main() -> None:
    bench_pairwise()


if __name__ == "__main__":
    main()
