"""
Demonstrates:
- Define precisely typed Python tasks using standard return type annotations.
- Retrieve strongly typed results using `Pipeline.run_typed` and a `TypedDict`
  schema that reflects channels kept in memory.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import TypedDict

from lahuta.pipeline import Pipeline, PipelineContext
from lahuta.sources import FileSource


# fmt: off
def file_info(ctx: PipelineContext) -> dict[str, str]:
    return {"file": os.path.basename(ctx.path), "ext": os.path.splitext(ctx.path)[1]}


def file_size(ctx: PipelineContext) -> dict[str, int]:
    return {"size": os.path.getsize(ctx.path)}


def display_name(ctx: PipelineContext) -> str:
    return f"processed: {os.path.basename(ctx.path)}"


class InfoRec(TypedDict):
    file: str
    ext:  str


class ResultsSchema(TypedDict):
    info:  list[InfoRec]
    sizes: list[dict[str, int]]
    names: list[str]


def main(files: list[str]) -> ResultsSchema:
    p = Pipeline(FileSource(files))

    p.add_task(name="info",  task=file_info)
    p.add_task(name="sizes", task=file_size)
    p.add_task(name="names", task=display_name)

    # IDEs and type checkers should know exact shapes per channel
    out = p.run_typed(ResultsSchema, threads=1)
    return out


if __name__ == "__main__":
    sample = [str(Path(__file__).resolve().parents[3] / "core" / "data" / "1kx2_small.cif")]
    results = main(sample)
    for rec in results["info"]:
        print(rec)
    for name in results["names"]:
        print(name)
