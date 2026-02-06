# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     class Email: pass
#     setattr(Email, "a", "besian")
#     setattr(Email, "b", "sejdiu")
#     setattr(Email, "c", "@gmail.com")
#     print(Email.a + Email.b + Email.c)
#
"""Entity bindings (EntityID/Kind) and records (AtomRec, RingRec, GroupRec), including RDKit Atom access."""

from __future__ import annotations

import pytest

from lahuta import EntityID, Kind, LahutaSystem, Topology, rdkit


@pytest.fixture(scope="session")
def topo(luni: LahutaSystem) -> Topology:
    assert luni.build_topology() is True
    return luni.get_topology()


def test_entities_and_records_rdkit(luni: LahutaSystem, topo: Topology) -> None:
    atom_recs = topo.atom_types
    assert isinstance(atom_recs, list) and len(atom_recs) == luni.n_atoms
    a0 = atom_recs[0]
    a0_idx = int(a0.idx())

    eid_atom = EntityID.make(Kind.Atom, a0_idx)
    assert eid_atom.kind == Kind.Atom and eid_atom.index == a0_idx
    assert str(eid_atom).startswith("Atom#")

    mol = topo.molecule()
    rd_atom = mol.getAtomWithIdx(a0_idx)
    assert isinstance(rd_atom, rdkit.Atom)
    assert rd_atom.getIdx() == a0_idx
    assert rd_atom.getAtomicNum() > 0

    rings = topo.rings
    assert isinstance(rings, list)
    if rings:
        r0 = rings[0]
        assert r0.size >= 3
        atoms = r0.atoms
        assert isinstance(atoms, list) and len(atoms) == r0.size
        assert all(hasattr(a, "getIdx") for a in atoms)

        eid_ring = EntityID.make(Kind.Ring, 0)
        assert str(eid_ring).startswith("Ring#0")

    groups = topo.groups
    assert isinstance(groups, list)
    if groups:
        g0 = groups[0]
        atoms = g0.atoms
        assert isinstance(atoms, list) and len(atoms) >= 1
        assert all(hasattr(a, "getIdx") for a in atoms)

        eid_group = EntityID.make(Kind.Group, 0)
        assert str(eid_group).startswith("Group#0")


def test_hypothesis_sampling_over_indices(luni: LahutaSystem, topo: Topology) -> None:
    hyp = pytest.importorskip("hypothesis", reason="Hypothesis not installed")
    st = pytest.importorskip("hypothesis.strategies", reason="Hypothesis not installed")

    atom_count = len(topo.atom_types)
    ring_count = len(topo.rings)
    group_count = len(topo.groups)

    @hyp.given(st.integers(min_value=0, max_value=max(0, atom_count - 1)))
    def _atoms(i: int) -> None:
        # RDKit Atom by index agrees with AtomRec.idx()
        rec = topo.get_atom(i)
        assert int(rec.idx()) == i
        mol = topo.molecule()
        a = mol.getAtomWithIdx(i)
        assert isinstance(a, rdkit.Atom) and a.getIdx() == i

    @hyp.given(st.integers(min_value=0, max_value=max(0, ring_count - 1)))
    def _rings(i: int) -> None:
        r = topo.get_ring(i)
        assert r.size >= 3
        atoms = r.atoms
        assert isinstance(atoms, list) and len(atoms) == r.size

    @hyp.given(st.integers(min_value=0, max_value=max(0, group_count - 1)))
    def _groups(i: int) -> None:
        g = topo.get_group(i)
        atoms = g.atoms
        assert isinstance(atoms, list) and len(atoms) >= 1

    _atoms()
    _rings()
    _groups()


def test_entity_id_kind_identity() -> None:
    """EntityID.kind returns the canonical Kind member (identity-safe)."""
    from lahuta import EntityID, Kind

    eid = EntityID.make(Kind.Ring, 42)
    assert int(eid) == (1 << 56) | 42
    assert eid.kind is Kind.Ring
    assert eid.kind == Kind.Ring
