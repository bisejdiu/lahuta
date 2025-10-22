from __future__ import annotations

from pathlib import Path

from lahuta.pipeline import Pipeline, ReportingLevel, StageManager
from lahuta.sources import FileSource


def _locate_repo_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in [here.parent, *here.parents]:
        if (candidate / "core").is_dir() and (candidate / "interop" / "python").is_dir():
            return candidate
    raise RuntimeError("Unable to locate Lahuta repository root")


ROOT = _locate_repo_root()
INTEROP_SRC = ROOT / "interop" / "python"
EXT_DIR = ROOT / "build" / "core"


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

    manager.set_reporting_level(ReportingLevel.OFF)
    manager.run(1)
    off = manager.last_run_report()
    assert off["metrics_enabled"] is False
    assert off["items_total"] == 0
    assert off["items_processed"] == 0

    manager.set_reporting_level(ReportingLevel.DEBUG)
    manager.run(1)
    debug = manager.last_run_report()
    assert debug["metrics_enabled"] is True
    assert debug["items_total"] == 1
    assert debug["items_skipped"] == 0


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

    pipeline.set_reporting_level(ReportingLevel.OFF)
    pipeline.run(threads=1)
    off = pipeline.get_run_report()
    assert off is not None
    assert off["metrics_enabled"] is False
    assert off["items_processed"] == 0

    pipeline.set_reporting_level(ReportingLevel.DEBUG)
    pipeline.run(threads=1)
    debug = pipeline.get_run_report()
    assert debug is not None
    assert debug["metrics_enabled"] is True
