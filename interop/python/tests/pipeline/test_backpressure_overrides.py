# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     class A: v = "besian"
#     class B(A): v = "sejdiu"
#     class C(B): v = "@gmail.com"
#     print("".join(c.v for c in reversed(C.__mro__[:-1])))
#
from __future__ import annotations

from lahuta.lib import lahuta as _lib
from lahuta.pipeline import BackpressureConfig, OnFull, get_default_backpressure_config, set_default_backpressure_config


def test_per_sink_backpressure_override_drops_large_payloads() -> None:
    # Default allows large payloads, and per-sink override enforces small queue budget.

    #
    # Note on flush timeout:
    #   StageManager closes the queue and then waits for writer threads in ChannelMultiplexer::close_and_flush.
    #   When a required sink rejects every payload (strict queue budget + DropLatest), the writer thread
    #   still has to observe close(), flush(), and report finished_. With the production default of 60s,
    #   tests would appear to "hang" for a minute even though they drop the payload immediately. I set
    #   the timeout to 5s here so the join guard does not dominate test runtimes while still exercising
    #   the backpressure logic.         - Besian, October 2025
    #
    original = get_default_backpressure_config()
    try:
        src = _lib.pipeline.sources.FileSource(["item"])
        sm = _lib.pipeline.StageManager(src)
        sm.set_flush_timeout(1.0)

        payload = "X" * 2048  # = 2 KiB

        def emit(_: _lib.pipeline.PipelineContext) -> str:
            return payload

        sm.add_python("emit", [], emit, "emit", serialize=True, store=True)

        ms = _lib.pipeline.MemorySink()
        cfg = BackpressureConfig()
        cfg.max_queue_bytes = 512
        cfg.max_batch_bytes = 512
        cfg.on_full = OnFull.DropLatest
        sm.connect_sink("emit", ms, cfg)

        sm.run(threads=1)
        assert ms.result() == []
    finally:
        set_default_backpressure_config(original)


def test_sink_reuse_keeps_original_config() -> None:
    # Reusing the same sink instance must keep the config captured at first connect
    original = get_default_backpressure_config()
    try:
        src = _lib.pipeline.sources.FileSource(["x"])
        sm = _lib.pipeline.StageManager(src)
        sm.set_flush_timeout(1.0)

        big = "Z" * 4096  # = 4 KiB

        def emit_a(_: _lib.pipeline.PipelineContext) -> str:
            return big

        def emit_b(_: _lib.pipeline.PipelineContext) -> str:
            return big

        sm.add_python("A", [], emit_a, "A", serialize=True, store=True)
        sm.add_python("B", [], emit_b, "B", serialize=True, store=True)

        ms = _lib.pipeline.MemorySink()

        # First connection with a strict config: small bytes, droplatest
        strict = BackpressureConfig()
        strict.max_queue_bytes = 512
        strict.max_batch_bytes = 512
        strict.on_full = OnFull.DropLatest
        sm.connect_sink("A", ms, strict)

        # Change defaults to permissive, then connect the SAME sink to another channel without cfg
        permissive = original
        permissive.max_queue_bytes = 8 * 1024 * 1024  # = 8 MiB
        permissive.max_batch_bytes = 8 * 1024 * 1024
        permissive.on_full = OnFull.Block
        set_default_backpressure_config(permissive)

        sm.connect_sink("B", ms)  # should reuse the existing ingress with 'strict' config

        sm.run(threads=1)
        assert ms.result() == []  # Both large payloads should have been dropped by the strict ingress
    finally:
        set_default_backpressure_config(original)
