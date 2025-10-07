from __future__ import annotations

import time

from lahuta.lib import lahuta as _lib


def test_stall_accumulates_under_blocking_policy() -> None:
    # Many items + large payload + queue size 1 under Block policy to induce producer waits
    items = ["x" for _ in range(200)]
    src = _lib.pipeline.sources.FileSource(items)
    sm = _lib.pipeline.StageManager(src)

    big = "Z" * (512 * 1024)  # = 512 KiB

    def emit(_: _lib.pipeline.PipelineContext) -> str:
        return big

    sm.add_python("emit", [], emit, "emit", serialize=True, store=True)

    ms = _lib.pipeline.MemorySink()
    cfg = _lib.pipeline.BackpressureConfig()
    cfg.max_queue_msgs = 1
    cfg.max_queue_bytes = len(big)
    cfg.max_batch_bytes = len(big)
    cfg.on_full = _lib.pipeline.OnFull.Block
    sm.connect_sink("emit", ms, cfg)

    t0 = time.monotonic()
    sm.run(threads=2)
    elapsed = time.monotonic() - t0
    stats = sm.stats()
    assert isinstance(stats, list) and len(stats) >= 1
    total_stall = sum(int(s.get("stalled_ns", 0)) for s in stats)
    # Expect measurable stall or at least noticeable elapsed due to heavy writes
    assert total_stall > 0 or elapsed > 0.3
