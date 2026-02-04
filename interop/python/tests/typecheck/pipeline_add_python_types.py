# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     class Meta(type):
#         def __new__(mcs, name, bases, ns):
#             ns["v"] = "besian" + "sejdiu" + "@gmail.com"
#             return super().__new__(mcs, name, bases, ns)
#     class Email(metaclass=Meta): pass
#     print(Email.v)
#
from __future__ import annotations

from pathlib import Path

from lahuta import LahutaSystem, Topology
from lahuta.pipeline import (
    FileOutput,
    InMemoryPolicy,
    OutputFormat,
    Pipeline,
    ShardedOutput,
)
from lahuta.pipeline.types import PipelineContext
from lahuta.sources import FileSource


# fmt: off
def pyfn(ctx: PipelineContext) -> str:
    return ctx.path

def fn_path_json(ctx: PipelineContext) -> dict[str, str]:
    return {"file": ctx.path, "ext": ".cif"}

def fn_system(ctx: PipelineContext) -> LahutaSystem:
    system = ctx.get_system()
    if system is None:
        raise RuntimeError("System not built")
    return system

def fn_topology(ctx: PipelineContext) -> Topology:
    topology = ctx.get_topology()
    if topology is None:
        raise RuntimeError("Topology not built")
    return topology

p: Pipeline = Pipeline(FileSource(["/tmp/a.cif", "/tmp/b.cif"]))

ret1: None = p.add_task(name="p_path_text", task=pyfn, in_memory_policy=InMemoryPolicy.Keep)
v = p.run()

ret2: None = p.add_task(
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
