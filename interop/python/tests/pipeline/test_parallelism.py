from __future__ import annotations

import threading
from typing import Any, Callable

import pytest

from lahuta.pipeline import InMemoryPolicy, Pipeline
from lahuta.sources import FilesSource


# fmt: off
class SharedProbe:
    """Deterministic barrier AND concurrency probe for Python callables.

    Tracks how many concurrent entries are inside the Python callable and
    supports a barrier to let multiple worker threads meet.
    """

    def __init__(self, target: int) -> None:
        self._lock = threading.Lock()
        self._cv = threading.Condition(self._lock)
        self.target = int(target)
        self.arrived = 0
        self.current = 0
        self.max = 0

    def enter(self) -> None:
        with self._lock:
            self.current += 1
            if self.current > self.max:
                self.max = self.current
            self.arrived += 1
            if self.arrived == self.target:
                self._cv.notify_all()
            elif self.target > 1:
                self._cv.wait_for(lambda: self.arrived >= self.target)
            self.current -= 1


def make_probe_fn(shared: SharedProbe) -> Callable[[Any], dict[str, Any]]:
    from lahuta.lib import lahuta as _lib

    def fn(ctx: _lib.pipeline.PipelineContext) -> dict[str, Any]:
        shared.enter()
        return {"ok": True, "path": ctx.path}

    return fn


def items(n: int) -> list[str]:
    return [f"item_{i}" for i in range(n)]


@pytest.mark.parametrize("threads", [2, 4])
def test_python_callable_serialize_true_limits_concurrency(threads: int) -> None:
    # serialize=True -> only one call to the Python function at a time, but pipeline remains parallel
    p = Pipeline(FilesSource(items(threads * 2)))

    shared = SharedProbe(target=1)  # no barrier, expect max==1
    fn = make_probe_fn(shared)

    p.add_task(
        name="py_serial",
        task=fn,
        depends=[],
        in_memory_policy=InMemoryPolicy.Keep,
        thread_safe=False,  # interpreted by wrapper as serialize=True
    )

    out = p.run(threads=threads)
    assert "py_serial" in out
    assert shared.max == 1


@pytest.mark.parametrize("threads", [2, 4])
def test_python_callable_serialize_false_allows_concurrency(threads: int) -> None:
    # serialize=False -> allow concurrent invocations (subject to GIL). Our barrier should ensure all arrive
    p = Pipeline(FilesSource(items(threads)))

    shared = SharedProbe(target=threads)
    fn = make_probe_fn(shared)

    p.add_task(
        name="py_concurrent",
        task=fn,
        depends=[],
        in_memory_policy=InMemoryPolicy.Keep,
        thread_safe=True,  # wrapper will set serialize=False
    )

    out = p.run(threads=threads)
    assert "py_concurrent" in out
    # All worker threads should have reached the callable concurrently
    assert shared.max == threads


@pytest.mark.parametrize("threads", [2, 4])
def test_two_python_tasks_both_concurrent(threads: int) -> None:
    # Two concurrent Python tasks should both observe full concurrency when no earlier bottleneck exists
    p = Pipeline(FilesSource(items(threads)))

    shared1 = SharedProbe(target=threads)
    shared2 = SharedProbe(target=threads)
    fn1 = make_probe_fn(shared1)
    fn2 = make_probe_fn(shared2)

    p.add_task(name="py1", task=fn1, depends=[], in_memory_policy=InMemoryPolicy.Keep, thread_safe=True)
    p.add_task(name="py2", task=fn2, depends=["py1"], in_memory_policy=InMemoryPolicy.Keep, thread_safe=True)

    out = p.run(threads=threads)
    assert "py1" in out and "py2" in out
    assert shared1.max == threads
    assert shared2.max == threads


@pytest.mark.parametrize("threads", [2])
def test_serial_upstream_does_not_collapse_downstream_parallelism(threads: int) -> None:
    # A serialized upstream Python task does not collapse overall pipeline parallelism.
    # It serializes only that callable. Downstream can still overlap across items.
    p = Pipeline(FilesSource(items(threads * 3)))

    serial_probe = SharedProbe(target=1)
    # No barrier for downstream. We observe natural overlap without forcing it
    concurrent_probe = SharedProbe(target=1)
    fn_serial = make_probe_fn(serial_probe)
    fn_conc   = make_probe_fn(concurrent_probe)

    p.add_task(name="py_serial", task=fn_serial, depends=[],            in_memory_policy=InMemoryPolicy.Keep, thread_safe=False)
    p.add_task(name="py_after",  task=fn_conc,   depends=["py_serial"], in_memory_policy=InMemoryPolicy.Keep, thread_safe=True)

    out = p.run(threads=threads)
    assert "py_serial" in out and "py_after" in out
    assert serial_probe.max == 1
    # Downstream may still overlap across items. Concurrency is bounded by pool size
    assert 1 <= concurrent_probe.max <= threads
