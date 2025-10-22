from __future__ import annotations

import tempfile
import time
from pathlib import Path

import pytest

from lahuta.pipeline import (
    BackpressureConfig,
    InMemoryPolicy,
    OutputFormat,
    Pipeline,
    PipelineContext,
)
from lahuta.pipeline.types import FileOutput, ShardedOutput
from lahuta.sources import FileSource


def test_multithreaded_writer_sinks() -> None:
    """Test that multi-threaded writer sinks work correctly with multiple writer threads."""
    items = [f"item_{i}" for i in range(50)]

    def slow_write_task(ctx: PipelineContext) -> str:
        # Simulate slow I/O operation
        time.sleep(0.01)
        return f"slow:{ctx.path}"

    def fast_write_task(ctx: PipelineContext) -> str:
        return f"fast:{ctx.path}"

    # Test with single writer thread (baseline)
    p_single = Pipeline(FileSource(items))
    p_single.add_task(
        name="slow_task",
        task=slow_write_task,
        in_memory_policy=InMemoryPolicy.Keep,
        writer_threads=1,
    )

    start_time = time.time()
    result_single = p_single.run(threads=4)
    single_duration = time.time() - start_time

    # Test with multiple writer threads
    p_multi = Pipeline(FileSource(items))
    p_multi.add_task(
        name="slow_task",
        task=slow_write_task,
        in_memory_policy=InMemoryPolicy.Keep,
        writer_threads=4,
    )

    start_time = time.time()
    result_multi = p_multi.run(threads=4)
    multi_duration = time.time() - start_time

    # Verify results contain the same items (order may differ due to threading)
    assert set(result_single["slow_task"]) == set(result_multi["slow_task"])
    assert len(result_single["slow_task"]) == len(items)
    assert len(result_multi["slow_task"]) == len(items)

    # Multi-threaded should generally be faster for slow I/O tasks
    # (though we don't enforce strict timing due to test environment variability)


def test_multithreaded_file_sinks() -> None:
    """Test multi-threaded writer with file sinks."""
    items = [f"item_{i}" for i in range(20)]

    def write_task(ctx: PipelineContext) -> dict:
        return {"path": ctx.path, "data": f"content_for_{ctx.path}"}

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Test with single writer thread
        p1 = Pipeline(FileSource(items))
        p1.add_task(
            name="task1",
            task=write_task,
            out=[],
            in_memory_policy=InMemoryPolicy.Keep,
            writer_threads=1,
        )

        output_file1 = Path(tmp_dir) / "output1.jsonl"
        p1.to_files("task1", path=output_file1, writer_threads=1)
        result1 = p1.run(threads=2)

        # Test with multiple writer threads
        p2 = Pipeline(FileSource(items))
        p2.add_task(
            name="task2",
            task=write_task,
            out=[],
            in_memory_policy=InMemoryPolicy.Keep,
            writer_threads=3,
        )

        output_file2 = Path(tmp_dir) / "output2.jsonl"
        p2.to_files("task2", path=output_file2, writer_threads=3)
        result2 = p2.run(threads=2)

        # Verify both in-memory and file outputs are correct
        assert len(result1["task1"]) == len(items)
        assert len(result2["task2"]) == len(items)

        # Both should have written identical content (up to task name differences)
        assert output_file1.exists()
        assert output_file2.exists()

        # Count lines in files
        lines1 = sum(1 for _ in output_file1.read_text().splitlines() if _.strip())
        lines2 = sum(1 for _ in output_file2.read_text().splitlines() if _.strip())
        assert lines1 == len(items)
        assert lines2 == len(items)


def test_mixed_writer_thread_configurations() -> None:
    """Test pipeline with different writer thread configurations per channel."""
    items = [f"item_{i}" for i in range(15)]

    def task_a(ctx: PipelineContext) -> str:
        return f"A:{ctx.path}"

    def task_b(ctx: PipelineContext) -> str:
        return f"B:{ctx.path}"

    with tempfile.TemporaryDirectory() as tmp_dir:
        p = Pipeline(FileSource(items))

        # Channel A: single writer thread
        p.add_task(
            name="task_a",
            task=task_a,
            in_memory_policy=InMemoryPolicy.Keep,
            writer_threads=1,
        )

        # Channel B: multiple writer threads
        p.add_task(
            name="task_b",
            task=task_b,
            in_memory_policy=InMemoryPolicy.Keep,
            writer_threads=3,
        )

        # Add file sinks with different writer thread counts
        output_a = Path(tmp_dir) / "output_a.jsonl"
        output_b = Path(tmp_dir) / "output_b.jsonl"

        p.to_files("task_a", path=output_a, writer_threads=1)
        p.to_files("task_b", path=output_b, writer_threads=3)

        result = p.run(threads=4)

        # Verify both channels have correct output
        assert len(result["task_a"]) == len(items)
        assert len(result["task_b"]) == len(items)

        # Verify file outputs
        assert output_a.exists()
        assert output_b.exists()

        lines_a = sum(1 for _ in output_a.read_text().splitlines() if _.strip())
        lines_b = sum(1 for _ in output_b.read_text().splitlines() if _.strip())
        assert lines_a == len(items)
        assert lines_b == len(items)


def test_backpressure_config_with_writer_threads() -> None:
    """Test that BackpressureConfig.writer_threads is exposed and works."""
    cfg = BackpressureConfig()
    assert cfg.writer_threads == 1  # default

    cfg.writer_threads = 4
    assert cfg.writer_threads == 4

    # Test validation
    cfg.validate()  # should not raise

    cfg.writer_threads = 0
    try:
        cfg.validate()
        assert False, "Expected validation error for writer_threads=0"
    except ValueError as e:
        assert "writer_threads" in str(e)


def test_writer_thread_validation_errors() -> None:
    """Ensure invalid writer thread counts raise helpful errors."""
    p = Pipeline(FileSource(["item"]))

    with pytest.raises(ValueError):
        p.add_task(
            name="invalid_task",
            task=lambda ctx: ctx.path,
            writer_threads=0,
        )

    with tempfile.TemporaryDirectory() as tmp_dir:
        with pytest.raises(TypeError):
            p.to_files("some_channel", path=Path(tmp_dir) / "out.jsonl", writer_threads=True)


def test_sharded_files_with_multiple_writers() -> None:
    """Test sharded file output with multiple writer threads."""
    items = [f"item_{i}" for i in range(25)]

    def write_task(ctx: PipelineContext) -> dict:
        return {"path": ctx.path, "timestamp": time.time()}

    with tempfile.TemporaryDirectory() as tmp_dir:
        out_dir = Path(tmp_dir) / "shards"
        sharded_out = ShardedOutput(out_dir=str(out_dir), shard_size=8)

        p = Pipeline(FileSource(items))
        p.add_task(
            name="sharded_task",
            task=write_task,
            out=[sharded_out],
            in_memory_policy=InMemoryPolicy.Keep,
            writer_threads=4,
        )

        result = p.run(threads=3)

        # Verify in-memory result
        assert len(result["sharded_task"]) == len(items)

        # Verify sharded files were created
        assert out_dir.exists()
        shard_files = list(out_dir.glob("*.ndjson"))  # Look for ndjson files
        assert len(shard_files) >= 3  # Should create at least 3 shards for 25 items with shard_size=8

        # Total lines across all shards
        total_lines = 0
        for shard_file in shard_files:
            lines = sum(1 for _ in shard_file.read_text().splitlines() if _.strip())
            total_lines += lines

        assert total_lines == len(items)
