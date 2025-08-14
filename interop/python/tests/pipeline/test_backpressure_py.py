from __future__ import annotations

import tempfile
import time
from pathlib import Path

from lahuta.pipeline import InMemoryPolicy, OutputFormat, Pipeline, PipelineContext


def _items(n: int) -> list[str]:
    return [f"item_{i}" for i in range(n)]


def test_memory_and_sharded_file_sinks_produce_expected_outputs() -> None:
    """Python-level sanity: one task, two sinks (memory + sharded files).

    Validates end-to-end that all N records are collected in memory and that
    a small shard_size produces multiple shards with the expected total lines.
    This exercises the Python wrapper -> StageManager -> async sinks path.
    """
    p = Pipeline.from_files(_items(7))

    def make_payload(ctx: PipelineContext):
        return f"rec:{ctx.path}"

    p.add_task(name="emit", task=make_payload, in_memory_policy=InMemoryPolicy.Keep)

    with tempfile.TemporaryDirectory() as tmp:
        out_dir = Path(tmp) / "shards"
        p.to_sharded_files("emit", out_dir=str(out_dir), fmt=OutputFormat.JSON, shard_size=3)

        out = p.run(threads=2)
        # memory sink received all
        assert "emit" in out and len(out["emit"]) == 7

        # sharded files produced (3,3,1 lines)
        files = p.file_outputs().get("emit", [])
        assert len(files) >= 3
        # Count total lines = 7
        total = 0
        for f in files:
            with open(f, "r", encoding="utf-8") as fh:
                total += sum(1 for _ in fh)
        assert total == 7


def test_large_payloads_do_not_crash_and_are_collected() -> None:
    """Stress the Python to C++ emission path with large strings.

    Ensures the wrapper and async sinks handle big payloads without deadlock or
    loss. This does not aim to saturate C++ byte queues (those knobs are not
    exposed in Python) - it's a safety/robustness check for the Python-facing API.
    """
    p = Pipeline.from_files(_items(20))

    big = "X" * (256 * 1024)  # 256 KiB per record

    def emit_big(ctx: PipelineContext):
        return big

    p.add_task(name="big", task=emit_big, in_memory_policy=InMemoryPolicy.Keep)
    out = p.run(threads=2)
    assert "big" in out and len(out["big"]) == 20
    # spot-check size of a few payloads
    assert len(out["big"][0]) == len(big)


def test_two_channels_with_file_and_memory_isolation() -> None:
    """Two independent channels: one spills to files, one stays in memory.

    Validates the Python wrapper wires two channels correctly and that having a
    file-backed sink on one channel doesn't impact the memory-only channel.
    """
    p = Pipeline.from_files(_items(9))

    def emit_a(ctx: PipelineContext):
        return f"A:{ctx.path}"

    def emit_b(ctx: PipelineContext):
        return f"B:{ctx.path}"

    # Channel A: memory + sharded files
    p.add_task(name="A", task=emit_a, in_memory_policy=InMemoryPolicy.Keep)
    with tempfile.TemporaryDirectory() as tmp:
        p.to_sharded_files("A", out_dir=str(Path(tmp) / "A"), fmt=OutputFormat.JSON, shard_size=4)

        # Channel B: memory only
        p.add_task(name="B", task=emit_b, in_memory_policy=InMemoryPolicy.Keep)

        out = p.run(threads=3)
        assert "A" in out and len(out["A"]) == 9
        assert "B" in out and len(out["B"]) == 9


def test_slow_file_sink_heavy_output_does_not_starve_memory_channel() -> None:
    """Heavy file output (large payloads) should not starve a memory-only channel.

    We emit large payloads to a sharded file sink (slow path) and small payloads
    to a memory sink (fast path) concurrently. We assert both channels collect
    all items and total runtime remains within a generous bound.
    """
    N = 60
    p = Pipeline.from_files(_items(N))

    small = "ok"
    large = "L" * (512 * 1024)  # 512 KiB per record

    def fast(ctx: PipelineContext):
        return f"F:{small}:{ctx.path}"

    def slow(ctx: PipelineContext):
        return large

    p.add_task(name="fast", task=fast, in_memory_policy=InMemoryPolicy.Keep)
    with tempfile.TemporaryDirectory() as tmp:
        p.to_sharded_files("fast", out_dir=str(Path(tmp) / "fast"), fmt=OutputFormat.JSON, shard_size=1000)
        p.add_task(name="slow", task=slow, in_memory_policy=InMemoryPolicy.Keep)
        p.to_sharded_files("slow", out_dir=str(Path(tmp) / "slow"), fmt=OutputFormat.JSON, shard_size=5)

        t0 = time.monotonic()
        out = p.run(threads=6)
        elapsed = time.monotonic() - t0

        assert "fast" in out and len(out["fast"]) == N
        assert "slow" in out and len(out["slow"]) == N
        # Not a perf test, just a starvation check under heavy IO.
        assert elapsed < 3.0


def test_slow_python_callable_in_second_task_allows_fast_first_task() -> None:
    """A slow Python task should not prevent a fast sibling from producing outputs.

    We register the fast task first (so it executes first per item) and a slow
    task second (sleep). We assert all outputs are present and that runtime
    stays under a generous bound to detect egregious serialization.
    """
    import time as _t

    N = 80
    p = Pipeline.from_files(_items(N))

    def fast(ctx: PipelineContext):
        return f"F:{ctx.path}"

    def slow(ctx: PipelineContext):
        _t.sleep(0.02)
        return f"S:{ctx.path}"

    # Fast first, then slow
    p.add_task(name="fast", task=fast, in_memory_policy=InMemoryPolicy.Keep, thread_safe=True)
    p.add_task(name="slow", task=slow, in_memory_policy=InMemoryPolicy.Keep, thread_safe=True)

    t0 = time.monotonic()
    out = p.run(threads=6)
    elapsed = time.monotonic() - t0

    assert "fast" in out and len(out["fast"]) == N
    assert "slow" in out and len(out["slow"]) == N
    # Aims to catch unintended global serialization
    assert elapsed < 2.0
