"""
Pipeline reporting example showcasing how to control metrics levels.

Run this module directly to execute a small pipeline over the bundled
example structures while toggling Lahuta's reporting levels.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from lahuta.pipeline import Pipeline, ReportingLevel
from lahuta.sources import FileSource


def _sample_files(limit: int = 2) -> list[str]:
    """Pick a couple of sample structures from the repository data directory."""
    data_dir = Path(__file__).resolve().parents[3] / "core" / "data"
    if not data_dir.exists():
        raise RuntimeError(f"Expected Lahuta data directory at {data_dir}")

    candidates: list[str] = []
    for pattern in ("*.cif", "*.cif.gz", "*.pdb"):
        matches = sorted(data_dir.glob(pattern))
        for path in matches:
            candidates.append(str(path))
            if len(candidates) >= limit:
                return candidates

    if not candidates:
        raise RuntimeError("No sample structures found in core/data")
    return candidates


def _baseline_task(ctx) -> str:
    """Store the base file name for each processed item."""
    name = Path(ctx.path).name
    ctx.set_text("basename", name)
    return name


def _run_with_level(pipe: Pipeline, level: ReportingLevel, threads: int = 2) -> dict:
    """Run the pipeline with the requested reporting level."""
    pipe.set_reporting_level(level)
    pipe.run(threads=threads)
    report = pipe.get_run_report()
    if report is None:
        raise RuntimeError("Pipeline.get_run_report returned None")
    return report


def _summarize(report: dict[str, object]) -> str:
    interesting = {
        "metrics_enabled": report["metrics_enabled"],
        "items_processed": report["items_processed"],
        "total_seconds": report["total_seconds"],
        "cpu_seconds": report["cpu_seconds"],
        "io_seconds": report["io_seconds"],
        "peak_inflight_items": report.get("peak_inflight_items"),
        "average_queue_depth": report.get("average_queue_depth"),
        "permit_waits": {
            "events": report.get("permit_wait_events"),
            "total_seconds": report.get("permit_wait_total_seconds"),
            "avg_seconds": report.get("permit_wait_avg_seconds"),
            "min_seconds": report.get("permit_wait_min_seconds"),
            "max_seconds": report.get("permit_wait_max_seconds"),
        },
        "mux": {
            "sinks": report.get("mux_sink_count"),
            "queue_peak": report.get("mux_queue_depth_peak"),
            "enqueued_msgs": report.get("mux_enqueued_msgs"),
            "written_msgs": report.get("mux_written_msgs"),
            "drops": report.get("mux_drops"),
        },
        "stage_breakdown": report.get("stage_breakdown"),
    }
    return json.dumps(interesting, indent=2)


def main(files: Iterable[str] | None = None) -> None:
    """Execute the example across all reporting levels."""
    selected = list(files) if files else _sample_files(5)
    print(f"Processing {len(selected)} items:")
    for path in selected:
        print(f"  - {path}")

    pipeline = Pipeline(FileSource(selected))
    pipeline.add_task(name="basename", task=_baseline_task, channel="basename", store=True)

    for level in (ReportingLevel.BASIC, ReportingLevel.OFF, ReportingLevel.DEBUG):
        report = _run_with_level(pipeline, level)
        print(f"\n[{level.name}] summary")
        print(_summarize(report))


if __name__ == "__main__":
    main()
