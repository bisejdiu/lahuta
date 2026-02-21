# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         def __init__(self, s=""): self.s = s
#         def __or__(self, other): return Email(self.s + other)
#         def __str__(self): return self.s
#     print(str(Email() | "besian" | "sejdiu" | "@gmail.com"))
#
"""Static typing exercises for the Pipeline wrapper API."""

from __future__ import annotations

import os
import sys
from typing import TypedDict

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from lahuta.pipeline import ErrorRecord, Pipeline, PipelineContext, PyTaskFn


def file_desc(ctx: PipelineContext) -> dict[str, str]:
    """A Python task returning a JSON object (dictionary)."""
    # Implementation irrelevant for type checking
    return {"file": ctx.path, "ext": ".cif"}


def file_name(ctx: PipelineContext) -> str:
    """A Python task returning text (string)."""
    return ctx.path


# Accept a callable with a precise return type (Protocol-based PyTaskFn)
def accepts_json_task(fn: PyTaskFn[dict[str, str]]) -> PyTaskFn[dict[str, str]]:
    return fn


_ = accepts_json_task(file_desc)


class DescRec(TypedDict):
    file: str
    ext: str


class ResultsSchema(TypedDict):
    desc: list[DescRec]
    names: list[str]
    boom: list[ErrorRecord]


def demonstrate_typed_results(p: Pipeline) -> None:
    """Demonstrate how users can strongly type the run() result via run_typed."""
    out = p.run_typed(ResultsSchema, threads=1)

    reveal_type(out)
    reveal_type(out["names"][0])
    reveal_type(out["desc"][0])
    reveal_type(out["desc"][0]["file"])
    reveal_type(out["boom"][0]["error"]["message"])
