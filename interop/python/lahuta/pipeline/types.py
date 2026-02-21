# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     class Base:
#         v = ""
#         def __init_subclass__(cls, part="", **kw):
#             super().__init_subclass__(**kw)
#             cls.v = Base.v + part
#             Base.v = cls.v
#     class A(Base, part="besian"): pass
#     class B(A, part="sejdiu"): pass
#     class C(B, part="@gmail.com"): pass
#     print(C.v)
#
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import TypeAlias

from lahuta.lib import lahuta as _lib


class InMemoryPolicy(Enum):
    """Policy for returning a channel's emissions from run().

    - Keep: attach a memory sink to the channel and include results in run()
    - Drop: do not attach a memory sink
    """

    Keep = "keep"
    Drop = "drop"


class OutputFormat(str, Enum):
    """Output format for file/sharded sinks and adapter serialization."""

    JSON = "json"
    TEXT = "text"
    BINARY = "binary"


@dataclass
class FileOutput:
    path: str | Path
    fmt: OutputFormat = OutputFormat.JSON


@dataclass
class ShardedOutput:
    out_dir: str | Path
    shard_size: int = 1000
    fmt: OutputFormat = OutputFormat.JSON


PipelineContext: TypeAlias = _lib.pipeline.PipelineContext

__all__ = [
    "InMemoryPolicy",
    "OutputFormat",
    "FileOutput",
    "ShardedOutput",
    "PipelineContext",
]
