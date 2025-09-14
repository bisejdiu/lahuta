from __future__ import annotations

from pathlib import Path

import pytest

import lahuta as lxx


# fmt: off
@pytest.fixture(scope="function")
def topo(ubi_cif: Path) -> lxx.Topology:
    sys = lxx.LahutaSystem(str(ubi_cif))
    assert sys.build_topology() is True
    return sys.get_topology()


@pytest.fixture(scope="function")
def contacts(topo: lxx.Topology):
    engine = lxx.MolStarContactsEngine()
    return engine.compute(topo)


def test_topology_resolve_returns_correct_record_types(topo: lxx.Topology, contacts) -> None:
    # Resolve the first 20 contacts via typed Topology resolvers and validate types/indices
    n = min(20, contacts.size())
    for i in range(n):
        c = contacts[i]
        # lhs
        if c.lhs.kind == lxx.Kind.Atom:
            rec = topo.resolve_atom(c.lhs)
            assert isinstance(rec, lxx.AtomRec)
            assert int(rec.idx()) == c.lhs.index
        elif c.lhs.kind == lxx.Kind.Ring:
            rec = topo.resolve_ring(c.lhs)
            assert isinstance(rec, lxx.RingRec)
            assert rec.size >= 1
        else:
            rec = topo.resolve_group(c.lhs)
            assert isinstance(rec, lxx.GroupRec)
            assert len(rec.atoms) >= 1
        # rhs
        if c.rhs.kind == lxx.Kind.Atom:
            rec = topo.resolve_atom(c.rhs)
            assert isinstance(rec, lxx.AtomRec)
            assert int(rec.idx()) == c.rhs.index
        elif c.rhs.kind == lxx.Kind.Ring:
            rec = topo.resolve_ring(c.rhs)
            assert isinstance(rec, lxx.RingRec)
            assert rec.size >= 1
        else:
            rec = topo.resolve_group(c.rhs)
            assert isinstance(rec, lxx.GroupRec)
            assert len(rec.atoms) >= 1


def test_entityresolver_resolve_contact_and_resolve_all(topo: lxx.Topology, contacts) -> None:
    resolver = lxx.EntityResolver(topo)

    c0 = contacts[0]
    lhs, rhs = resolver.resolve_contact(c0)
    if c0.lhs.kind == lxx.Kind.Atom:
        assert isinstance(lhs, lxx.AtomRec) and int(lhs.idx()) == c0.lhs.index
    elif c0.lhs.kind == lxx.Kind.Ring:
        assert isinstance(lhs, lxx.RingRec) and lhs.size >= 1
    else:
        assert isinstance(lhs, lxx.GroupRec) and len(lhs.atoms) >= 1

    if c0.rhs.kind == lxx.Kind.Atom:
        assert isinstance(rhs, lxx.AtomRec) and int(rhs.idx()) == c0.rhs.index
    elif c0.rhs.kind == lxx.Kind.Ring:
        assert isinstance(rhs, lxx.RingRec) and rhs.size >= 1
    else:
        assert isinstance(rhs, lxx.GroupRec) and len(rhs.atoms) >= 1

    pairs = resolver.resolve_all(contacts)
    assert isinstance(pairs, list) and len(pairs) == contacts.size()


def test_non_atom_contact_subset_and_first_five_atom_pairs(topo: lxx.Topology, contacts) -> None:
    resolver = lxx.EntityResolver(topo)

    atom_atom = []
    non_atom_pairs = []
    for c in contacts:
        lhs, rhs = resolver.resolve_contact(c)
        lhs_is_atom = isinstance(lhs, lxx.AtomRec)
        rhs_is_atom = isinstance(rhs, lxx.AtomRec)
        if lhs_is_atom and rhs_is_atom:
            if len(atom_atom) < 5:
                atom_atom.append((int(lhs.idx()), int(rhs.idx()), c.type.category))
        else:
            non_atom_pairs.append((lhs, rhs, c))

    # First 5 atom-atom contacts: (lhs_idx, rhs_idx, Category)
    expected_first5 = [
        (4, 76,   lxx.Category.HydrogenBond),
        (8, 72,   lxx.Category.HydrogenBond),
        (23, 581, lxx.Category.HydrogenBond),
        (24, 61,  lxx.Category.HydrogenBond),
        (28, 37,  lxx.Category.HydrogenBond),
    ]
    assert atom_atom == expected_first5

    # Non-atom-involving contacts: 4 group-group Ionic contacts with expected sizes and charge polarity
    assert len(non_atom_pairs) == 4
    expected_sizes = [(3, 2), (2, 2), (1, 2), (1, 2)]
    expected_d2 = [23.597, 24.836, 23.156, 10.084]
    for i, (lhs, rhs, c) in enumerate(non_atom_pairs):
        # kinds
        assert isinstance(lhs, lxx.GroupRec) and isinstance(rhs, lxx.GroupRec)
        # interaction category
        assert c.type.category == lxx.Category.Ionic
        # sizes
        la, lb = expected_sizes[i]
        assert len(lhs.atoms) == la and len(rhs.atoms) == lb
        # charge polarity (Positive on lhs, Negative on rhs)
        assert lhs.a_type.has(lxx.AtomType.PositiveCharge)
        assert rhs.a_type.has(lxx.AtomType.NegativeCharge)
        # distance squared approx
        assert float(c.distance_sq) == pytest.approx(expected_d2[i], rel=0, abs=0.02)
