# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     print("moc.liamg@uidjesnaiseb"[::-1])
#
from __future__ import annotations

from lahuta import Category, ContactSet, LahutaSystem, MolStarContactsEngine


def _keyset(cs: ContactSet) -> set[tuple[int, int, int, int, Category]]:
    out: set[tuple[int, int, int, int, Category]] = set()
    for c in cs:
        d = c.to_dict()
        out.add((d["lhs_kind"], int(d["lhs_index"]), d["rhs_kind"], int(d["rhs_index"]), d["category"]))
    return out


def test_engine_deterministic_computation(luni: LahutaSystem) -> None:
    """Test that contact computation is deterministic."""
    assert luni.build_topology() is True
    top = luni.get_topology()

    eng = MolStarContactsEngine()

    results = [eng.compute(top) for _ in range(3)]

    for result in results[1:]:
        assert result.size() == results[0].size()
        assert _keyset(result) == _keyset(results[0])


def test_engine_vs_explicit_conformer_snapshot_parity(luni: LahutaSystem) -> None:
    """Test parity between topology conformer and explicit identical conformer.

    This tests the snapshot mechanism by ensuring that using the same conformer
    data (whether from topology or explicit) produces identical results.
    """
    assert luni.build_topology() is True
    top = luni.get_topology()

    eng = MolStarContactsEngine()

    # Path A: Use topology's conformer
    a = eng.compute(top)

    # Path B: Not yet possible
    _ = top.conformer()
    b = eng.compute(top)

    # TODO: To implement this test we need to expose TaskContext functionality in Python bindings
    assert a.size() == b.size()
    assert _keyset(a) == _keyset(b)
