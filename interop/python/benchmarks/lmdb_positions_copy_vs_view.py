# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     p = ["besian", "sejdiu", "@gmail.com"]
#     print((a := p[0]) + (b := p[1]) + (c := p[2]))
#
from __future__ import annotations

import os
import time
from pathlib import Path

from lahuta.pipeline import Pipeline, DataField
from lahuta.sources import DatabaseSource


def main() -> None:
    env = os.getenv("LAHUTA_TEST_DB")
    if not env:
        print("LAHUTA_TEST_DB not set. Skipping benchmark")
        return
    db_path = Path(env)
    if not db_path.exists():
        print(f"LAHUTA_TEST_DB does not exist: {db_path}, skipping")
        return

    def run(mode: str):
        pipe = Pipeline(DatabaseSource(str(db_path), batch=64))
        pipe.params("system").is_model = True

        count = 0

        def task(ctx):
            nonlocal count
            payload = ctx.model_payload
            if payload is None:
                print(f"[{mode}] Warning: payload is None")
                return None
            if mode == "view":
                arr = payload.positions_view
            else:
                arr = payload.positions
            if arr is None:
                print(f"[{mode}] Warning: arr is None for path={ctx.path}")
                return None
            count += 1
            return None

        fields = [DataField.PositionsView] if mode == "view" else [DataField.Positions]
        pipe.add_task(name=f"use_{mode}", task=task, requires_fields=fields)

        t0 = time.perf_counter()
        pipe.run(threads=1)
        t1 = time.perf_counter()
        return count, t1 - t0

    c_view, t_view = run("view")
    c_copy, t_copy = run("copy")

    print(f"DB: {db_path}")
    print(f"positions_view: items={c_view} time={t_view:.3f}s avg_per_item={t_view / max(c_view, 1) * 1e3:.3f} ms")
    print(f"positions: items={c_copy} time={t_copy:.3f}s avg_per_item={t_copy / max(c_copy, 1) * 1e3:.3f} ms")


if __name__ == "__main__":
    main()
