from __future__ import annotations

from pathlib import Path

from lahuta.pipeline import InMemoryPolicy, Pipeline, PipelineContext
from lahuta.sources import FilesSource


# fmt: off
def test_context_path_and_emit_text(data_path) -> None:
    data = str(data_path("ubi.cif"))

    def basename(ctx: PipelineContext) -> str:
        return Path(ctx.path).name

    p = Pipeline(FilesSource([data]))
    p.add_task(name="names", task=basename, depends=["system"], in_memory_policy=InMemoryPolicy.Keep)
    out = p.run(threads=1)
    assert out["names"] == [Path(data).name]


def test_context_system_topology_and_emit_json(data_path) -> None:
    data = str(data_path("ubi.cif"))

    def check(ctx: PipelineContext) -> dict:
        sys = ctx.get_system()
        top = ctx.get_topology()
        return {"have_sys": sys is not None, "have_top": top is not None}

    p = Pipeline(FilesSource([data]))
    # Depend on 'topology' to make sure a topology is present before running
    p.add_task(name="chk", task=check, depends=["topology"], in_memory_policy=InMemoryPolicy.Keep)
    out = p.run(threads=1)
    assert out["chk"] == [{"have_sys": True, "have_top": True}]


def test_context_read_upstream_text_and_skip_on_missing(data_path) -> None:
    """Regression for store vs context-read semantics.

    - With store=False on task 'a', its payload is not written to TaskContext.
      Downstream 'b' sees ctx.get_text('a') is None and returns None to skip emission.
    - With store=True on 'a', the payload is written under key 'a'. Downstream
      'b' reads it and emits.

    This overlaps with pipeline/test_pipeline_dynamic_store_vs_emit.py, which
    exercises the same contract in a more focused way. Keeping this here ties
    the behavior to the PipelineContext API examples and guards against API-level
    regressions (e.g., get_text None behavior, emission skipping on None).
    """
    data = str(data_path("ubi.cif"))

    def a(_: PipelineContext) -> str:
        return "A"

    def b(ctx: PipelineContext):
        s = ctx.get_text("a")
        if s is None:
            return None
        return f"B:{s}"

    p = Pipeline(FilesSource([data]))
    p.add_task(name="a", task=a, depends=["system"], in_memory_policy=InMemoryPolicy.Keep, store=False)
    p.add_task(name="b", task=b, depends=["a"],      in_memory_policy=InMemoryPolicy.Keep)
    out = p.run(threads=1)
    assert out["a"] == ["A"]
    assert out.get("b", []) == []

    # Now enable store to allow B to read A's payload
    p = Pipeline(FilesSource([data]))
    p.add_task(name="a", task=a, depends=["system"], in_memory_policy=InMemoryPolicy.Keep, store=True)
    p.add_task(name="b", task=b, depends=["a"],      in_memory_policy=InMemoryPolicy.Keep)
    out = p.run(threads=1)
    assert out["a"] == ["A"]
    assert out["b"] == ["B:A"]


def test_return_none_skips_emission(data_path) -> None:
    data = str(data_path("ubi.cif"))

    def noop(_: PipelineContext):
        return None

    p = Pipeline(FilesSource([data]))
    p.add_task(name="n", task=noop)
    out = p.run(threads=1)
    assert out.get("n", []) == []


def test_context_frame_metadata_defaults(data_path) -> None:
    data = str(data_path("ubi.cif"))

    captured: dict[str, object] = {}

    def capture(ctx: PipelineContext) -> str:
        captured["conformer_id"] = ctx.conformer_id
        captured["session_id"] = ctx.session_id
        captured["timestamp_ps"] = ctx.timestamp_ps
        captured["meta"] = ctx.frame_metadata()
        return "ok"

    p = Pipeline(FilesSource([data]))
    p.add_task(
        name="capture",
        task=capture,
        depends=["system"],
        in_memory_policy=InMemoryPolicy.Keep,
    )
    out = p.run(threads=1)

    assert out["capture"] == ["ok"]
    assert captured["conformer_id"] == 0
    assert captured["session_id"] == data
    assert captured["timestamp_ps"] is None

    meta = captured["meta"]
    assert isinstance(meta, dict)
    assert meta["session_id"] == data
    assert meta["conformer_id"] == 0
    assert meta["timestamp_ps"] is None
