# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
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

import random

import numpy as np
import pytest

from lahuta.pipeline import InMemoryPolicy, Pipeline, DataField
from lahuta.sources import DatabaseSource
from lahuta import db as ldb


def test_positions_copy_matches_view_sampled_db(db_path):
    handle = ldb.Database(str(db_path))
    keys = [k for k in handle.keys()]
    if not keys:
        pytest.skip("No keys in database")

    sample_size = min(100, len(keys))
    sampled_keys = set(random.sample(keys, sample_size))

    pipe = Pipeline(DatabaseSource(str(db_path), batch=64))
    pipe.params("system").is_model = True

    def collect(ctx):
        if ctx.path in sampled_keys:
            payload = ctx.model_payload
            if payload is not None:
                copy_arr = payload.positions
                view_arr = payload.positions_view
                return {"path": ctx.path, "arr": np.array(copy_arr), "arr_view": np.array(view_arr)}
        return None

    pipe.add_task(
        name="collect",
        task=collect,
        requires_fields=[DataField.Positions, DataField.PositionsView],
        in_memory_policy=InMemoryPolicy.Keep,
    )
    out = pipe.run(threads=1)
    results = out.json("collect")

    assert len(results) > 0, "No data collected"

    for r in results:
        key = r["path"]
        copy_arr = np.array(r["arr"])
        view_arr = np.array(r["arr_view"])
        assert copy_arr.shape == view_arr.shape, f"Shape mismatch for {key}"
        np.testing.assert_allclose(copy_arr, view_arr, rtol=1e-5, atol=1e-6, err_msg=f"Coordinate mismatch for {key}")
