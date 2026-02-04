# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     p = ["besian", "sejdiu", "@gmail.com"]
#     print((a := p[0]) + (b := p[1]) + (c := p[2]))
#
"""Exercise configurable flush timeout behaviour exposed to Python."""

from __future__ import annotations

import time
from typing import Iterable

import pytest

from lahuta.lib import lahuta as _lib


def _make_stage_manager(channel: str, items: Iterable[str]) -> _lib.pipeline.StageManager:
    """Create a StageManager emitting payloads to the requested channel."""
    src = _lib.pipeline.sources.FileSource(list(items))
    sm = _lib.pipeline.StageManager(src)
    sm.set_auto_builtins(False)

    def emit(ctx: _lib.pipeline.PipelineContext) -> str:
        return f"payload:{ctx.path}"

    sm.add_python("emit", [], emit, channel, serialize=True, store=False)
    return sm


def test_required_sink_timeout_raises() -> None:
    sm = _make_stage_manager("slow", [f"item{i}" for i in range(2)])
    slow_sink = _lib.pipeline.SlowSink(sleep_seconds=0.0, flush_sleep_seconds=0.2)
    cfg = _lib.pipeline.BackpressureConfig()
    cfg.required = True
    sm.connect_sink("slow", slow_sink, cfg)
    sm.set_flush_timeout(0.01)

    with pytest.raises(RuntimeError, match="required sink failed: sink deadline exceeded"):
        sm.run(threads=1)


def test_non_required_sink_timeout_does_not_raise() -> None:
    sm = _make_stage_manager("optional", ["only"])
    slow_sink = _lib.pipeline.SlowSink(sleep_seconds=0.0, flush_sleep_seconds=0.05)
    cfg = _lib.pipeline.BackpressureConfig()
    cfg.required = False
    sm.connect_sink("optional", slow_sink, cfg)
    sm.set_flush_timeout(0.001)

    sm.run(threads=1)
    time.sleep(0.06)
    assert slow_sink.count() == 1


def test_flush_timeout_roundtrip_and_validation() -> None:
    sm = _make_stage_manager("unused", ["item"])
    # Default should match production value
    assert sm.get_flush_timeout() == pytest.approx(60.0)

    # Tests run faster if we lower the timeout explicitly
    sm.set_flush_timeout(5.0)
    assert sm.get_flush_timeout() == pytest.approx(5.0)

    sm.set_flush_timeout(1.5)
    assert sm.get_flush_timeout() == pytest.approx(1.5)

    with pytest.raises(ValueError):
        sm.set_flush_timeout(-0.1)


def test_run_flush_timeout_kwarg_restores_previous() -> None:
    sm = _make_stage_manager("slow_kwarg", [f"item{i}" for i in range(2)])
    slow_sink = _lib.pipeline.SlowSink(sleep_seconds=0.0, flush_sleep_seconds=0.2)
    cfg = _lib.pipeline.BackpressureConfig()
    cfg.required = True
    sm.connect_sink("slow_kwarg", slow_sink, cfg)
    sm.set_flush_timeout(0.5)

    with pytest.raises(RuntimeError, match="required sink failed: sink deadline exceeded"):
        sm.run(threads=1, flush_timeout=0.01)

    assert sm.get_flush_timeout() == pytest.approx(0.5)


def test_run_flush_timeout_kwarg_validation() -> None:
    sm = _make_stage_manager("noop", ["item"])
    with pytest.raises(ValueError):
        sm.run(threads=1, flush_timeout=-1.0)
