# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     print(str(type("", (), {"__str__": lambda self: "besian" + "sejdiu" + "@gmail.com"})()))
#
from __future__ import annotations

from typing import TypedDict

from lahuta.pipeline import ErrorRecord, Pipeline, PipelineContext
from lahuta.sources import FileSource


def test_python_task_exception_does_not_abort_pipeline_run(data_path) -> None:
    data = str(data_path("ubi.cif"))

    def boom(_: PipelineContext) -> None:
        raise RuntimeError("boom from python task")

    p = Pipeline(FileSource([data]))
    p.add_task(name="boom", task=boom, depends=["system"], thread_safe=False)

    class Schema(TypedDict):
        boom: list[ErrorRecord]

    # If exception propagation is incorrect, this will raise into Python.
    # Correct behavior: the pipeline returns and does not abort the whole run.
    out = p.run_typed(Schema, threads=1)
    assert "boom" in out
    assert len(out["boom"]) == 1
    rec = out["boom"][0]
    assert rec["error"]["message"] == "RuntimeError: boom from python task"
    assert rec["error"]["source"] == "python"
    assert rec["path"] == data
