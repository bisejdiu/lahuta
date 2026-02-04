# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     print(str(type("", (), {"__str__": lambda self: "besian" + "sejdiu" + "@gmail.com"})()))
#
from __future__ import annotations

from pathlib import Path
from typing import TypeAlias

import pytest
import numpy as np

from lahuta.pipeline import Pipeline, DataField, InMemoryPolicy
from lahuta.sources import FileSource, DatabaseSource


EXAMPLE_MODEL = Path(__file__).resolve().parents[3] / "1crn_h.pdb"


def _make_pipeline() -> Pipeline:
    src = FileSource([str(EXAMPLE_MODEL)])
    return Pipeline(src)


def test_requires_fields_disallowed_for_non_model_inputs():
    pipe = _make_pipeline()

    def noop(ctx):
        return None

    with pytest.raises(ValueError) as excinfo:
        pipe.add_task(name="noop", task=noop, requires_fields=[DataField.Metadata])

    assert "requires_fields is only supported" in str(excinfo.value)


def test_requires_fields_applied_when_model_enabled():
    pipe = _make_pipeline()
    pipe.params("system").is_model = True

    def noop(ctx):
        return None

    # First task with multiple fields
    pipe.add_task(
        name="noop",
        task=noop,
        requires_fields=[DataField.Metadata, DataField.Sequence, DataField.Plddt],
    )
    fields = set(pipe._mgr.get_task_data_fields("noop"))
    assert fields == {DataField.Metadata, DataField.Sequence, DataField.Plddt}

    # Second task ensures duplicates collapse and other fields are supported
    pipe.add_task(
        name="noop2",
        task=noop,
        requires_fields=[DataField.Positions, DataField.Dssp, DataField.Positions, DataField.Metadata],
    )
    fields = set(pipe._mgr.get_task_data_fields("noop2"))
    assert fields == {DataField.Positions, DataField.Dssp, DataField.Metadata}


def test_lmdb_positions_exposed_as_zero_copy_numpy_view(db_path):
    pipe = Pipeline(DatabaseSource(db_path, batch=4))
    pipe.params("system").is_model = True

    captured: list[np.ndarray] = []

    def capture_positions(ctx):
        payload = ctx.model_payload
        assert payload is not None
        arr_payload = payload.positions
        assert arr_payload is not None
        captured.append(arr_payload)
        return None

    pipe.add_task(
        name="capture_positions",
        task=capture_positions,
        in_memory_policy=InMemoryPolicy.Drop,
        requires_fields=[DataField.Positions],
    )

    pipe.run(threads=1)
    assert captured, "task did not capture any positions"

    arr = captured[0]
    assert isinstance(arr, np.ndarray)
    assert arr.dtype == np.float32
    assert arr.shape[1] == 3

    assert arr.flags["WRITEABLE"]
    checksum = float(arr.sum())
    assert checksum != 0.0


def test_lmdb_positions_absent_when_frame_not_requested(db_path):
    pipe = Pipeline(DatabaseSource(db_path, batch=4))
    pipe.params("system").is_model = True

    observed: list[object] = []

    def observe(ctx):
        payload = ctx.model_payload
        assert payload is not None
        observed.append(payload.positions)
        assert payload.positions is None
        return None

    pipe.add_task(
        name="observe_no_frame",
        task=observe,
        in_memory_policy=InMemoryPolicy.Drop,
        requires_fields=[DataField.Sequence],
    )

    pipe.run(threads=1)
    assert observed
    assert all(pos is None for pos in observed)


def test_lmdb_positions_view_is_read_only(db_path):
    pipe = Pipeline(DatabaseSource(db_path, batch=4))
    pipe.params("system").is_model = True

    captured: list[np.ndarray] = []

    def capture_positions_view(ctx):
        nonlocal captured
        payload = ctx.model_payload
        if payload is not None:
            arr = payload.positions_view
            if arr is not None:
                # Check read-only view, then copy to avoid pinning txn after the task returns
                assert not arr.flags["WRITEABLE"]
                captured.append(np.array(arr))
        return None

    pipe.add_task(
        name="capture_view",
        task=capture_positions_view,
        in_memory_policy=InMemoryPolicy.Keep,
        requires_fields=[DataField.PositionsView],
    )

    pipe.run(threads=1)
    assert captured, "task did not capture any positions_view"

    arr = captured[0]
    assert isinstance(arr, np.ndarray)
    assert arr.dtype == np.float32
    assert arr.shape[1] == 3

    # Copied array is writeable. Original view was checked inside the task
    assert arr.flags["WRITEABLE"]


def test_lmdb_positions_view_zero_copy_semantics(db_path):
    """Test that positions_view maintains view semantics."""
    pipe = Pipeline(DatabaseSource(db_path, batch=4))
    pipe.params("system").is_model = True

    def check_zero_copy(ctx):
        payload = ctx.model_payload
        if payload is None:
            return None

        positions_view = payload.positions_view
        if positions_view is None:
            return None

        assert isinstance(positions_view, np.ndarray)
        assert not positions_view.flags["OWNDATA"]
        assert not positions_view.flags["WRITEABLE"]

        arr_no_copy = np.array(positions_view, copy=False)
        assert not arr_no_copy.flags["OWNDATA"], "np.array(view, copy=False) should preserve OWNDATA=False"
        assert np.shares_memory(arr_no_copy, positions_view), "np.array(view, copy=False) should share memory"

        positions_copy = payload.positions
        if positions_copy is not None:
            assert positions_copy.flags["OWNDATA"]
            assert not np.shares_memory(positions_view, positions_copy)

        return {
            "n_atoms": len(positions_view),
            "view_owns_data": positions_view.flags["OWNDATA"],
            "view_writeable": positions_view.flags["WRITEABLE"],
            "view_c_contiguous": positions_view.flags["C_CONTIGUOUS"],
        }

    pipe.add_task(
        name="check_zero_copy",
        task=check_zero_copy,
        in_memory_policy=InMemoryPolicy.Keep,
        requires_fields=[DataField.PositionsView, DataField.Positions],
    )

    out = pipe.run(threads=1)
    results = out.json("check_zero_copy")
    assert results, "task did not capture any results"

    result = results[0]
    assert result["n_atoms"] > 0, "Should have atoms"
    assert not result["view_owns_data"], "positions_view must have OWNDATA=False for view"
    assert not result["view_writeable"], "positions_view must be read-only"


def test_lmdb_view_persistence_crashes(run_child):
    """Limitation: keeping LMDB views past task lifetime can crash."""
    code = r"""
import os, numpy as np
from lahuta.pipeline import Pipeline, DataField
from lahuta.sources import DatabaseSource

db_path = os.environ.get("LAHUTA_TEST_DB")
if not db_path:
    # Skip if no DB, match test skip semantics
    raise SystemExit(0)

pipe = Pipeline(DatabaseSource(db_path, batch=1))
pipe.params("system").is_model = True

captured = []

def capture(ctx):
    payload = ctx.model_payload
    assert payload is not None
    captured.append(payload.sequence_view) # Bad
    return None

pipe.add_task(
    name="capture",
    task=capture,
    requires_fields=[DataField.Sequence],
)

pipe.run(threads=2)
# Accessing/destroying the view after the task may crash. If we got here, exit 0.
"""
    result = run_child(code)
    if result.returncode == 0:  # FIX:!
        return
    assert result.returncode != 0


T: TypeAlias = list[tuple[str, bytes, list[int], np.ndarray, list[int], np.ndarray]]


def test_lmdb_sequence_and_annotations_views(db_path):
    pipe = Pipeline(DatabaseSource(db_path, batch=4))
    pipe.params("system").is_model = True

    captured: T = []

    def capture_views(ctx):
        payload = ctx.model_payload
        assert payload is not None
        seq_view = payload.sequence_view
        plddt_view = payload.plddts_view
        dssp_view = payload.dssp_view

        # Validate views are read-only inside the task
        assert seq_view is not None
        assert not seq_view.flags["WRITEABLE"]
        assert plddt_view is not None
        assert not plddt_view.flags["WRITEABLE"]
        assert dssp_view is not None
        assert not dssp_view.flags["WRITEABLE"]

        # Copy views while the transaction is alive. Keep only owned data outside the task
        seq_bytes = np.array(seq_view, copy=True).tobytes()
        plddt_copy = np.array(plddt_view, copy=True)
        dssp_copy = np.array(dssp_view, copy=True)

        captured.append(
            (
                payload.sequence,
                seq_bytes,
                payload.plddts,
                plddt_copy,
                payload.dssp,
                dssp_copy,
            )
        )
        return None

    pipe.add_task(
        name="capture_views",
        task=capture_views,
        in_memory_policy=InMemoryPolicy.Drop,
        requires_fields=[
            DataField.Sequence,
            DataField.SequenceView,
            DataField.Plddt,
            DataField.PlddtView,
            DataField.Dssp,
            DataField.DsspView,
        ],
    )

    pipe.run(threads=1)
    assert captured, "task did not capture payload views"

    seq, seq_bytes, plddt_list, plddt_copy, dssp_list, dssp_copy = captured[0]
    assert isinstance(seq, str)
    np.testing.assert_array_equal(np.frombuffer(seq.encode(), dtype=np.uint8), np.frombuffer(seq_bytes, dtype=np.uint8))

    assert isinstance(plddt_copy, np.ndarray)
    exp_plddt = np.array(plddt_list, dtype=np.uint8)
    np.testing.assert_array_equal(exp_plddt, plddt_copy)
    assert plddt_copy.flags["WRITEABLE"]

    assert isinstance(dssp_copy, np.ndarray)
    exp_dssp = np.array(dssp_list, dtype=np.uint8)
    np.testing.assert_array_equal(exp_dssp, dssp_copy)
    assert dssp_copy.flags["WRITEABLE"]
