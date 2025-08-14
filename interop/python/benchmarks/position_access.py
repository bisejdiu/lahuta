from __future__ import annotations

import gc
import math
import os
import statistics as stats
import time
from pathlib import Path
from timeit import Timer
from typing import Callable, Iterable

from lahuta.lib import lahuta as lxx

# fmt: off

def _ns_stats(samples_ns: Iterable[float]) -> dict[str, float]:
    lst = list(samples_ns); lst.sort()
    n   = len(lst)

    mean  = stats.fmean(lst)  if n else float("nan")
    stdev = stats.pstdev(lst) if n else float("nan")
    med = lst[n // 2] if n else float("nan")
    p95 = lst[math.floor(0.95 * (n - 1))] if n else float("nan")
    p99 = lst[math.floor(0.99 * (n - 1))] if n else float("nan")
    return {
        "mean":   mean,
        "stdev":  stdev,
        "min":    lst[0] if n else float("nan"),
        "median": med,
        "p95":    p95,
        "p99":    p99,
        "max":    lst[-1] if n else float("nan"),
    }


def benchmark_callable(name: str, fn: Callable[[], object], *, number: int = 1000, repeat: int = 15, baseline_call_ns: float | None = None) -> dict[str, float]:
    """Benchmark a zero-arg callable using timeit, returning per-call ns stats.

    - Uses timeit.Timer with repeat, number to amortize timer overhead.
    - GC is disabled by timeit internally. We keep it off around warmup too.
    - Optionally substracts a provided baseline per-call overhead (e.g., a no-op).
    """
    # Warmup outside measurement to pre-import modules and trigger caches
    was_enabled = gc.isenabled()
    if was_enabled:
        gc.disable()
    try:
        for _ in range(5):
            fn()
    finally:
        if was_enabled:
            gc.enable()

    t = Timer(fn)
    # timeit disables GC during timing. Returns total seconds for 'number' calls.
    totals_s = t.repeat(repeat=repeat, number=number)
    per_call_ns = [(total / number) * 1e9 for total in totals_s]

    if baseline_call_ns is not None:
        per_call_ns = [max(0.0, x - baseline_call_ns) for x in per_call_ns]

    res = _ns_stats(per_call_ns)
    print(f"\n{name}  (per-call, ns)  [{number=}, {repeat=}]")
    print(
        f"  mean={res['mean']:.1f}  stdev={res['stdev']:.1f}  median={res['median']:.1f}  "
        f"min={res['min']:.1f}  p95={res['p95']:.1f}  p99={res['p99']:.1f}  max={res['max']:.1f}"
    )
    return res


def main() -> None:
    DATA   = Path(__file__).resolve().parents[3] / "core"/ "data" / "fubi.cif"
    number = int(os.environ.get("LAHUTA_BENCH_NUMBER", "1000"))
    repeat = int(os.environ.get("LAHUTA_BENCH_REPEAT", "20"))
    sys    = lxx.LahutaSystem.from_model_file(str(DATA))
    props  = sys.props

    # first-call latencies before any warmups
    def _first_call_ns(fn: Callable[[], object]) -> float:
        was_enabled = gc.isenabled()
        gc.collect()
        if was_enabled:
            gc.disable()
        try:
            t0 = time.perf_counter_ns()
            _ = fn()
            t1 = time.perf_counter_ns()
        finally:
            if was_enabled:
                gc.enable()
        return float(t1 - t0)

    # Here order DOES matter! finicky, finicky...
    cold_pv = _first_call_ns(lambda: props.positions_view)
    cold_pc = _first_call_ns(lambda: props.positions)
    print("\nCold first call (one-shot) latencies:")
    print(f"  props.positions_view: {cold_pv:.1f} ns")
    print(f"  props.positions     : {cold_pc:.1f} ns")

    # Baseline: empty Python callable overhead per call
    base_ns = benchmark_callable("baseline(lambda: None)", lambda: None, number=number, repeat=repeat)["mean"]

    # Baseline: lightweight property call (int) to approximate descriptor+bind overhead
    int_prop = benchmark_callable(
        "sys.n_atoms (int property baseline)",
        lambda: sys.n_atoms,
        number=number,
        repeat=repeat,
        baseline_call_ns=base_ns,
    )
    int_prop_mean = int_prop["mean"]

    # Measure positions_view
    pv = benchmark_callable(
        "props.positions_view",
        lambda: props.positions_view,
        number=number,
        repeat=repeat,
        baseline_call_ns=base_ns,
    )
    print(f"  Delta vs int-prop baseline: {pv['mean'] - int_prop_mean:.1f} ns (mean)")

    # Measure positions
    pc = benchmark_callable(
        "props.positions",
        lambda: props.positions,
        number=number,
        repeat=repeat,
        baseline_call_ns=base_ns,
    )
    print(f"  Delta vs int-prop baseline: {pc['mean'] - int_prop_mean:.1f} ns (mean)")


if __name__ == "__main__":
    main()
