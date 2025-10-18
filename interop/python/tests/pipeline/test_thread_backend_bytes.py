from __future__ import annotations

from pathlib import Path

import pytest

from lahuta.pipeline import InMemoryPolicy, Pipeline, PipelineContext
from lahuta.sources import FileSource

DATA_PATH = Path("core/data/models/AF-P0CL56-F1-model_v4.cif.gz")


def emit_bytes(_: PipelineContext) -> bytes:
    return b"\x00\x01hello\xff"


def write_blob(ctx: PipelineContext) -> None:
    ctx.set_bytes("blob", b"\xaa\xbb")


def echo_blob(ctx: PipelineContext) -> bytes:
    data = ctx.get_bytes("blob")
    assert data is not None
    return data + b"\xcc"


@pytest.mark.skipif(not DATA_PATH.exists(), reason="Test data file not found")
def test_thread_backend_binary_payload_channel() -> None:
    p = Pipeline(FileSource([str(DATA_PATH)]))
    p.add_task(name="bin", task=emit_bytes, in_memory_policy=InMemoryPolicy.Keep, store=False)

    res = p.run(threads=1, backend="threads")
    raw = res.raw("bin")
    assert isinstance(raw, tuple) and len(raw) == 1
    assert raw[0] == b"\x00\x01hello\xff"

    with pytest.raises(ValueError):
        res.strings("bin")


@pytest.mark.skipif(not DATA_PATH.exists(), reason="ubi.cif fixture not available")
def test_thread_backend_store_bytes_roundtrip() -> None:
    p = Pipeline(FileSource([str(DATA_PATH)]))
    p.add_task(name="write_blob", task=write_blob, in_memory_policy=InMemoryPolicy.Drop, store=False)
    p.add_task(
        name="echo_blob", task=echo_blob, depends=["write_blob"], in_memory_policy=InMemoryPolicy.Keep, store=False
    )

    res = p.run(threads=1, backend="threads")
    raw = res.raw("echo_blob")
    assert raw == (b"\xaa\xbb\xcc",)
