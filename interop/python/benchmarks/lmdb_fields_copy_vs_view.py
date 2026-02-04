# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     print(functools.reduce(lambda a, b: a + b, ["besian", "sejdiu", "@gmail.com"], ""))
#
from __future__ import annotations

import os
import time
from pathlib import Path

from lahuta.pipeline import Pipeline, DataField, ReportingLevel
from lahuta.sources import DatabaseSource


def _bench(db_path: Path, label: str, fields: list[DataField], accessor_name: str) -> tuple[int, float]:
    pipe = Pipeline(DatabaseSource(str(db_path), batch=64))
    pipe.params("system").is_model = True
    pipe.set_reporting_level(ReportingLevel.OFF)

    count = 0

    def task(ctx):
        nonlocal count
        payload = ctx.model_payload
        if payload is None:
            return None
        arr = getattr(payload, accessor_name)
        if arr is None:
            return None
        count += 1
        return None

    pipe.add_task(name=f"use_{label}", task=task, requires_fields=fields)
    t0 = time.perf_counter()
    pipe.run(threads=1)
    t1 = time.perf_counter()
    return count, t1 - t0


def main() -> None:
    env = os.getenv("LAHUTA_TEST_DB")
    if not env:
        print("LAHUTA_TEST_DB not set. Skipping benchmark")
        return
    db_path = Path(env)
    if not db_path.exists():
        print(f"LAHUTA_TEST_DB does not exist: {db_path}, skipping")
        return

    results: list[tuple[str, int, float]] = []

    # Sequence
    c, t = _bench(db_path, "sequence_view", [DataField.SequenceView], "sequence_view")
    results.append(("sequence_view", c, t))
    c, t = _bench(db_path, "sequence_copy", [DataField.Sequence], "sequence")
    results.append(("sequence_copy", c, t))

    # pLDDT
    c, t = _bench(db_path, "plddt_view", [DataField.PlddtView], "plddts_view")
    results.append(("plddt_view", c, t))
    c, t = _bench(db_path, "plddt_copy", [DataField.Plddt], "plddts")
    results.append(("plddt_copy", c, t))

    # DSSP
    c, t = _bench(db_path, "dssp_view", [DataField.DsspView], "dssp_view")
    results.append(("dssp_view", c, t))
    c, t = _bench(db_path, "dssp_copy", [DataField.Dssp], "dssp")
    results.append(("dssp_copy", c, t))

    print(f"DB: {db_path}")
    for name, c, t in results:
        per = t / max(c, 1) * 1e3
        print(f"{name:>14}: items={c:6d} time={t:8.3f}s avg_per_item={per:8.3f} ms")


if __name__ == "__main__":
    main()
