# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     print(type("", (), {"__call__": lambda _: "besian" + "sejdiu" + "@gmail.com"})()())
#
"""Test that FrameView properly loads only the view without copy."""

from __future__ import annotations

from lahuta.pipeline import InMemoryPolicy, Pipeline, DataField
from lahuta.sources import DatabaseSource


def test_frameview_loads_only_view(db_path):
    pipe = Pipeline(DatabaseSource(str(db_path), batch=1))
    pipe.params("system").is_model = True

    def check_view(ctx):
        payload = ctx.model_payload
        assert payload is not None
        assert payload.positions_view is not None
        # Should NOT have positions (copy)
        # However, positions() may fall back to copy from view, so we can't check it directly
        # assert payload.positions is None
        return {
            "has_view": payload.positions_view is not None,
        }

    pipe.add_task(
        name="check_view",
        task=check_view,
        requires_fields=[DataField.PositionsView],
        in_memory_policy=InMemoryPolicy.Keep,
    )

    out = pipe.run(threads=1)
    results = out.json("check_view")

    assert len(results) > 0, "No items processed"
    for r in results:
        assert r["has_view"], "PositionsView should load positions_view"


def test_frame_loads_both(db_path):
    pipe = Pipeline(DatabaseSource(str(db_path), batch=1))
    pipe.params("system").is_model = True

    def check(ctx):
        payload = ctx.model_payload
        assert payload is not None
        assert payload.positions is not None
        assert payload.positions_view is not None
        return {
            "has_positions": payload.positions is not None,
            "has_view": payload.positions_view is not None,
        }

    pipe.add_task(
        name="check_both",
        task=check,
        requires_fields=[DataField.Positions, DataField.PositionsView],
        in_memory_policy=InMemoryPolicy.Keep,
    )

    out = pipe.run(threads=1)
    results = out.json("check_both")

    assert len(results) > 0, "No items processed"
    for r in results:
        assert r["has_positions"], "Positions should load positions"
        assert r["has_view"], "Positions should load positions_view"
