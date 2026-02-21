# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         def __repr__(self):
#             return "besian" + "sejdiu" + "@gmail.com"
#     print(repr(Email()))
#
from __future__ import annotations

from pathlib import Path

import pytest

from lahuta import (
    AtomRec,
    AtomType,
    Category,
    EntityResolver,
    GroupRec,
    Kind,
    LahutaSystem,
    MolStarContactsEngine,
    RingRec,
    Topology,
)


# fmt: off
@pytest.fixture(scope="function")
def topo(ubi_cif: Path) -> Topology:
    sys = LahutaSystem(str(ubi_cif))
    assert sys.build_topology() is True
    return sys.get_topology()


@pytest.fixture(scope="function")
def contacts(topo: Topology):
    engine = MolStarContactsEngine()
    return engine.compute(topo)


def test_topology_resolve_returns_correct_record_types(topo: Topology, contacts) -> None:
    # Resolve the first 20 contacts via typed Topology resolvers and validate types/indices
    n = min(20, contacts.size())
    for i in range(n):
        c = contacts[i]
        # lhs
        if c.lhs.kind == Kind.Atom:
            rec = topo.resolve_atom(c.lhs)
            assert isinstance(rec, AtomRec)
            assert int(rec.idx) == c.lhs.index
        elif c.lhs.kind == Kind.Ring:
            rec = topo.resolve_ring(c.lhs)
            assert isinstance(rec, RingRec)
            assert rec.size >= 1
        else:
            rec = topo.resolve_group(c.lhs)
            assert isinstance(rec, GroupRec)
            assert len(rec.atoms) >= 1
        # rhs
        if c.rhs.kind == Kind.Atom:
            rec = topo.resolve_atom(c.rhs)
            assert isinstance(rec, AtomRec)
            assert int(rec.idx) == c.rhs.index
        elif c.rhs.kind == Kind.Ring:
            rec = topo.resolve_ring(c.rhs)
            assert isinstance(rec, RingRec)
            assert rec.size >= 1
        else:
            rec = topo.resolve_group(c.rhs)
            assert isinstance(rec, GroupRec)
            assert len(rec.atoms) >= 1


def test_entityresolver_resolve_contact_and_resolve_all(topo: Topology, contacts) -> None:
    resolver = EntityResolver(topo)

    c0 = contacts[0]
    lhs, rhs = resolver.resolve_contact(c0)
    if c0.lhs.kind == Kind.Atom:
        assert isinstance(lhs, AtomRec) and int(lhs.idx) == c0.lhs.index
    elif c0.lhs.kind == Kind.Ring:
        assert isinstance(lhs, RingRec) and lhs.size >= 1
    else:
        assert isinstance(lhs, GroupRec) and len(lhs.atoms) >= 1

    if c0.rhs.kind == Kind.Atom:
        assert isinstance(rhs, AtomRec) and int(rhs.idx) == c0.rhs.index
    elif c0.rhs.kind == Kind.Ring:
        assert isinstance(rhs, RingRec) and rhs.size >= 1
    else:
        assert isinstance(rhs, GroupRec) and len(rhs.atoms) >= 1

    pairs = resolver.resolve_all(contacts)
    assert isinstance(pairs, list) and len(pairs) == contacts.size()


def test_non_atom_contact_subset_and_first_five_atom_pairs(topo: Topology, contacts) -> None:
    resolver = EntityResolver(topo)

    atom_atom = []
    non_atom_pairs = []
    for c in contacts:
        lhs, rhs = resolver.resolve_contact(c)
        lhs_is_atom = isinstance(lhs, AtomRec)
        rhs_is_atom = isinstance(rhs, AtomRec)
        if lhs_is_atom and rhs_is_atom:
            if len(atom_atom) < 5:
                atom_atom.append((int(lhs.idx), int(rhs.idx), c.type.category))
        else:
            non_atom_pairs.append((lhs, rhs, c))

    # First 5 atom-atom contacts: (lhs_idx, rhs_idx, Category)
    expected_first5 = [
        (4, 76,   Category.HydrogenBond),
        (8, 72,   Category.HydrogenBond),
        (23, 581, Category.HydrogenBond),
        (24, 61,  Category.HydrogenBond),
        (28, 37,  Category.HydrogenBond),
    ]
    assert atom_atom == expected_first5

    # Non-atom-involving contacts: 4 group-group Ionic contacts with expected sizes and charge polarity
    assert len(non_atom_pairs) == 4
    expected_sizes = [(3, 2), (2, 2), (1, 2), (1, 2)]
    expected_d2 = [23.597, 24.836, 23.156, 10.084]
    for i, (lhs, rhs, c) in enumerate(non_atom_pairs):
        # kinds
        assert isinstance(lhs, GroupRec) and isinstance(rhs, GroupRec)
        # interaction category
        assert c.type.category == Category.Ionic
        # sizes
        la, lb = expected_sizes[i]
        assert len(lhs.atoms) == la and len(rhs.atoms) == lb
        # charge polarity (Positive on lhs, Negative on rhs)
        assert lhs.a_type.has(AtomType.PositiveCharge)
        assert rhs.a_type.has(AtomType.NegativeCharge)
        # distance squared approx
        assert float(c.distance_sq) == pytest.approx(expected_d2[i], rel=0, abs=0.02)
