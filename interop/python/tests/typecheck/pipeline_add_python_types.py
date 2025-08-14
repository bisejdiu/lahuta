from __future__ import annotations

from pathlib import Path

from lahuta.lib import lahuta as lxx
from lahuta.pipeline import (
    FileOutput,
    InMemoryPolicy,
    OutputFormat,
    Pipeline,
    ShardedOutput,
)
from lahuta.pipeline.types import PipelineContext


# fmt: off
def pyfn(ctx: PipelineContext) -> str:
    return ctx.path

def fn_path_json(ctx: PipelineContext) -> dict[str, str]:
    return {"file": ctx.path, "ext": ".cif"}

def fn_system(ctx: PipelineContext) -> lxx.LahutaSystem:
    system = ctx.get_system()
    if system is None:
        raise RuntimeError("System not built")
    return system

def fn_topology(ctx: PipelineContext) -> lxx.Topology:
    topology = ctx.get_topology()
    if topology is None:
        raise RuntimeError("Topology not built")
    return topology

p: Pipeline = Pipeline.from_files(["/tmp/a.cif", "/tmp/b.cif"])

_p1: Pipeline = p.add_task(name="p_path_text", task=pyfn, in_memory_policy=InMemoryPolicy.Keep)
v = _p1.run()

_p2: Pipeline = p.add_task(
    name="kwargs_all",
    task=pyfn,
    channel="chan",
    thread_safe=False,
    store=True,
    depends=["system", "topology"],
    out=[
        FileOutput(Path("file.ndjson"), fmt=OutputFormat.JSON),
        ShardedOutput(Path("shards"), shard_size=10, fmt=OutputFormat.JSON),
    ],
)
