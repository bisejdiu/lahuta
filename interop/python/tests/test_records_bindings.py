"""
Verifies:
- AtomRec  indices cover all atoms and are within bounds.
- RingRec  invariants (size >= 3, atom indices valid, center/normal numeric).
- GroupRec invariants (atoms non-empty, indices valid, center numeric).
- Property access and basic mutability for exposed fields.
"""

from __future__ import annotations

import numpy as np
import pytest

import lahuta as lxx

from pathlib import Path


# fmt: off
@pytest.fixture(scope="function")
def topo(ubi_cif: Path) -> lxx.Topology:
    # Build to Complete to expose all records deterministically
    sys = lxx.LahutaSystem(str(ubi_cif))
    assert sys.build_topology() is True
    return sys.get_topology()


def test_atomrec_indices_cover_atoms(luni: lxx.LahutaSystem, topo: lxx.Topology) -> None:
    recs = topo.atom_types
    n = luni.n_atoms
    assert len(recs) == n
    idxs = [int(r.idx()) for r in recs]
    assert min(idxs) == 0 and max(idxs) == n - 1
    # It's possible order is 0..n-1, check uniqueness
    assert len(set(idxs)) == n


def _is_point3d(obj) -> bool:
    return all(hasattr(obj, f) for f in ("x", "y", "z"))


def test_ringrec_invariants(luni: lxx.LahutaSystem, topo: lxx.Topology) -> None:
    n = luni.n_atoms
    rings = topo.rings

    assert isinstance(rings, list)
    assert len(rings) >= 1
    for r in rings:
        assert r.size >= 3
        atoms = r.atoms
        assert isinstance(atoms, list) and len(atoms) == r.size
        for a in atoms:
            idx = int(a.getIdx())
            assert 0 <= idx < n
        assert isinstance(r.aromatic, (bool, np.bool_))
        assert _is_point3d(r.center) and _is_point3d(r.normal)


def test_grouprec_invariants(luni: lxx.LahutaSystem, topo: lxx.Topology) -> None:
    n = luni.n_atoms
    groups = topo.groups
    assert isinstance(groups, list)
    assert len(groups) >= 1
    for g in groups:
        atoms = g.atoms
        assert isinstance(atoms, list)
        assert len(atoms) >= 1
        for a in atoms:
            idx = int(a.getIdx())
            assert 0 <= idx < n
        assert _is_point3d(g.center)


def test_atomrec_property_setter_no_crash(topo: lxx.Topology) -> None:
    # Setting to the same value should be a no-op and not crash.
    recs = topo.atom_types
    if not recs:
        pytest.skip("No atom records exposed")
    rec0 = recs[0]
    t0 = int(rec0.type)
    rec0.type = rec0.type  # set to same value
    assert int(rec0.type) == t0


def test_hypothesis_indexed_access(topo: lxx.Topology) -> None:
    hyp = pytest.importorskip("hypothesis",            reason="Hypothesis not installed")
    st  = pytest.importorskip("hypothesis.strategies", reason="Hypothesis not installed")

    atom_count  = len(topo.atom_types)
    ring_count  = len(topo.rings)
    group_count = len(topo.groups)

    @hyp.given(st.integers(min_value=0, max_value=max(0, atom_count - 1)))
    def _atom_idx_prop(i: int) -> None:
        r = topo.get_atom(i)
        assert int(r.idx()) == i

    @hyp.given(st.integers(min_value=0, max_value=max(0, ring_count - 1)))
    def _ring_idx_prop(i: int) -> None:
        r = topo.get_ring(i)
        assert r.size >= 3

    @hyp.given(st.integers(min_value=0, max_value=max(0, group_count - 1)))
    def _group_idx_prop(i: int) -> None:
        g = topo.get_group(i)
        assert len(g.atoms) >= 1

    _atom_idx_prop()
    _ring_idx_prop()
    _group_idx_prop()
