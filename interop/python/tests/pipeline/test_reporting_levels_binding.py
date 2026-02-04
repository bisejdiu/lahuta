# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     async def f():
#         return "besian" + "sejdiu" + "@gmail.com"
#     print(asyncio.run(f()))
#
from __future__ import annotations

from pathlib import Path

from lahuta.pipeline import Pipeline, ReportingLevel, StageManager
from lahuta.sources import FileSource


def test_stage_manager_reporting_levels(data_path):
    sample = str(data_path("ubi.cif"))
    manager = StageManager(FileSource([sample]))
    manager.set_auto_builtins(False)

    def record(ctx):
        return Path(ctx.path).name

    manager.add_python("basename", fn=record, serialize=False, store=False)

    manager.set_reporting_level(ReportingLevel.BASIC)
    manager.run(1)
    basic = manager.last_run_report()
    assert basic["metrics_enabled"] is True
    assert basic["items_total"] == 1
    assert basic["items_processed"] == 1
    assert basic["peak_inflight_items"] >= 1
    assert basic["permit_wait_events"] == 0
    assert basic["stage_breakdown"] == []
    assert "mux_sink_count" in basic

    manager.set_reporting_level(ReportingLevel.OFF)
    manager.run(1)
    off = manager.last_run_report()
    assert off["metrics_enabled"] is False
    assert off["items_total"] == 0
    assert off["items_processed"] == 0
    assert off["peak_inflight_items"] == 0
    assert off["permit_wait_events"] == 0
    assert off["stage_breakdown"] == []
    assert "mux_sink_count" in off

    manager.set_reporting_level(ReportingLevel.DEBUG)
    manager.run(1)
    debug = manager.last_run_report()
    assert debug["metrics_enabled"] is True
    assert debug["items_total"] == 1
    assert debug["items_skipped"] == 0
    assert debug["peak_inflight_items"] >= 1
    assert debug["permit_wait_events"] == 0
    assert isinstance(debug["stage_breakdown"], list)
    assert len(debug["stage_breakdown"]) == debug["stage_count"]
    if debug["stage_breakdown"]:
        first_stage = debug["stage_breakdown"][0]
        assert "label" in first_stage and "compute_seconds" in first_stage
    assert "mux_sink_count" in debug


def test_pipeline_reporting_levels(data_path):
    sample = str(data_path("ubi.cif"))
    pipeline = Pipeline(FileSource([sample]))

    def record(ctx):
        return Path(ctx.path).name

    pipeline.add_task(name="basename", task=record, store=False)

    pipeline.set_reporting_level(ReportingLevel.BASIC)
    pipeline.run(threads=1)
    basic = pipeline.get_run_report()
    assert basic is not None
    assert basic["metrics_enabled"] is True
    assert basic["items_processed"] == 1
    assert basic["peak_inflight_items"] >= 1
    assert basic["stage_breakdown"] == []
    assert "mux_sink_count" in basic

    pipeline.set_reporting_level(ReportingLevel.OFF)
    pipeline.run(threads=1)
    off = pipeline.get_run_report()
    assert off is not None
    assert off["metrics_enabled"] is False
    assert off["items_processed"] == 0
    assert off["peak_inflight_items"] == 0
    assert "mux_sink_count" in off

    pipeline.set_reporting_level(ReportingLevel.DEBUG)
    pipeline.run(threads=1)
    debug = pipeline.get_run_report()
    assert debug is not None
    assert debug["metrics_enabled"] is True
    assert isinstance(debug["stage_breakdown"], list)
    assert len(debug["stage_breakdown"]) == debug["stage_count"]
    assert "mux_sink_count" in debug
