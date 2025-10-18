from __future__ import annotations

from pathlib import Path

from lahuta.pipeline import InMemoryPolicy, Pipeline, PipelineContext
from lahuta.sources import FileSource


# fmt: off
def _stage_one(ctx: PipelineContext) -> dict[str, int | str]:
    topology = ctx.get_topology()
    assert topology is not None

    atom_ids = list(topology.get_atom_ids())
    atom_count = len(atom_ids)

    ctx.set_text("shared_atom_count", str(atom_count))
    return {"atoms": atom_count, "path": ctx.path}


def _stage_two(ctx: PipelineContext) -> str:
    value = ctx.get_text("shared_atom_count")
    assert value is not None
    return value


def test_process_backend_preserves_context_and_matches_threads(data_path) -> None:
    """Ensure the process backend executes via StageManager and mirrors threaded results."""
    data_file = Path(data_path("ubi.cif"))
    assert data_file.exists(), "Test dataset 'ubi.cif' must be available"

    pipeline = Pipeline(FileSource([str(data_file)]))

    pipeline.add_task(
        name="stage_one",
        task=_stage_one,
        depends=["topology"],
        in_memory_policy=InMemoryPolicy.Keep,
    )
    pipeline.add_task(
        name="stage_two",
        task=_stage_two,
        depends=["stage_one"],
        in_memory_policy=InMemoryPolicy.Keep,
    )

    threaded_first  = pipeline.run(threads=1, backend="threads")
    process_result  = pipeline.run(threads=1, backend="processes", processes=2)
    threaded_second = pipeline.run(threads=1, backend="threads")

    assert threaded_first .get("stage_two") == process_result.get("stage_two")
    assert threaded_second.get("stage_two") == process_result.get("stage_two")
    assert threaded_first .get("stage_one")["atoms"] == int(process_result.get("stage_two"))


def test_process_backend_falls_back_for_local_callables(data_path) -> None:
    """Local/nested callables should execute correctly by falling back to the threaded backend."""
    data_file = Path(data_path("ubi.cif"))
    assert data_file.exists(), "Test dataset 'ubi.cif' must be available"

    def local_stage(_: PipelineContext) -> str:
        return "local"

    pipeline = Pipeline(FileSource([str(data_file)]))
    pipeline.add_task(
        name="local",
        task=local_stage,
        depends=["topology"],
        in_memory_policy=InMemoryPolicy.Keep,
    )

    threaded  = pipeline.run(threads=1, backend="threads")
    processes = pipeline.run(threads=1, backend="processes")

    assert threaded.get("local") == processes.get("local")
