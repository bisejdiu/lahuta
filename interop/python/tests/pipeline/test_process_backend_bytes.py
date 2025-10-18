from __future__ import annotations

from pathlib import Path

import pytest

from lahuta.pipeline import InMemoryPolicy, Pipeline, PipelineContext
from lahuta.sources import FileSource

# fmt: off
DATA_PATH = Path("core/data/models/AF-P0CL56-F1-model_v4.cif.gz")


def emit_bytes(_: PipelineContext) -> bytes:
    return b"\x00\x01hello\xff"


def write_blob(ctx: PipelineContext) -> None:
    # Set a binary blob in the worker context store
    ctx.set_bytes("blob", b"\x10\x11\x12")


def echo_blob(ctx: PipelineContext) -> bytes:
    data = ctx.get_bytes("blob")
    assert data is not None
    return data + b"!"


@pytest.mark.skipif(not DATA_PATH.exists(), reason="Test data file not found")
def test_process_backend_binary_payload_channel() -> None:
    p = Pipeline(FileSource([str(DATA_PATH)]))
    p.add_task(name="bin", task=emit_bytes, in_memory_policy=InMemoryPolicy.Keep, store=False)

    res = p.run(threads=1, backend="processes", processes=1)
    raw = res.raw("bin")
    assert isinstance(raw, tuple) and len(raw) == 1
    assert raw[0] == b"\x00\x01hello\xff"

    with pytest.raises(ValueError):
        res.strings("bin")


@pytest.mark.skipif(not DATA_PATH.exists(), reason="Test data file not found")
def test_process_backend_store_bytes_roundtrip() -> None:
    p = Pipeline(FileSource([str(DATA_PATH)]))
    p.add_task(name="write_blob", task=write_blob, in_memory_policy=InMemoryPolicy.Drop, store=False)
    p.add_task(name="echo_blob", task=echo_blob, depends=["write_blob"], in_memory_policy=InMemoryPolicy.Keep, store=False)

    res = p.run(threads=1, backend="processes", processes=1)
    raw = res.raw("echo_blob")
    assert raw == (b"\x10\x11\x12!",)
