# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     print(functools.reduce(lambda a, b: a + b, ["besian", "sejdiu", "@gmail.com"], ""))
#
"""Tests for entity-based contact detection."""

from __future__ import annotations

from typing import Iterable

import pytest

from lahuta import (
    AtomType,
    Category,
    Contact,
    InteractionType,
    InteractionTypeSet,
    Kind,
    LahutaSystem,
    Topology,
)
from lahuta.entities import atoms, compute_contacts, find_contacts, groups, rings


# fmt: off
@pytest.fixture(scope="session")
def topo(luni: LahutaSystem) -> Topology:
    assert luni.build_topology() is True
    return luni.get_topology()


def _key(c: Contact) -> tuple[Kind, int, Kind, int, Category]:
    d = c.to_dict()
    return d["lhs_kind"], int(d["lhs_index"]), d["rhs_kind"], int(d["rhs_index"]), d["category"]


def test_system_sanity(luni: LahutaSystem, topo: Topology) -> None:
    assert isinstance(luni.n_atoms, int) and luni.n_atoms == 671
    assert len(topo.atom_records) == 671
    assert len(topo.rings)  >= 1
    assert len(topo.groups) >= 1


def test_provider_engine_molstar_counts_and_samples(luni: LahutaSystem, topo: Topology) -> None:
    ms_all = compute_contacts(topo, provider="molstar")
    hbonds = compute_contacts(topo, provider="molstar", only=InteractionType.HydrogenBond)

    assert ms_all.size() == 157
    assert hbonds.size() == 102

    # Spot-check exemplars via indices and category using structured data
    keys_all = {_key(c) for c in ms_all}
    keys_hb  = {_key(c) for c in hbonds}

    def atom_pair(i: int, j: int, cat: Category) -> tuple[Kind, int, Kind, int, Category]:
        return (Kind.Atom, min(i, j), Kind.Atom, max(i, j), cat)

    assert atom_pair(4,  76, Category.HydrogenBond) in keys_all
    assert atom_pair(8,  72, Category.HydrogenBond) in keys_all
    assert atom_pair(28, 52, Category.WeakHydrogenBond) in keys_all
    assert atom_pair(29, 37, Category.HydrogenBond) in keys_all

    assert atom_pair(4, 76,  Category.HydrogenBond) in keys_hb
    assert atom_pair(8, 72,  Category.HydrogenBond) in keys_hb
    assert atom_pair(29, 37, Category.HydrogenBond) in keys_hb


def test_compute_contacts_multiple_interactions(topo: Topology) -> None:
    hb  = InteractionType.HydrogenBond
    whb = InteractionType.WeakHydrogenBond

    combo = compute_contacts(topo, provider="molstar", only=[hb, whb])
    assert combo.size() > 0

    combo_contacts = set(combo)
    categories = {c.type.category for c in combo_contacts}
    assert Category.HydrogenBond in categories
    assert Category.WeakHydrogenBond in categories
    assert categories.issubset({Category.HydrogenBond, Category.WeakHydrogenBond})

    combo_or = compute_contacts(topo, provider="molstar", only=hb | whb)
    assert set(combo_or) == combo_contacts

    filter_set = InteractionTypeSet(hb)
    filter_set |= whb
    combo_set = compute_contacts(topo, provider="molstar", only=filter_set)
    assert set(combo_set) == combo_contacts

    nested = compute_contacts(topo, provider="molstar", only=[hb, [whb, hb], (InteractionTypeSet(whb),)])
    assert set(nested) == combo_contacts

    def pair_generator() -> Iterable[InteractionType]:
        yield from (hb, whb)

    combo_iter = compute_contacts(topo, provider="molstar", only=pair_generator())
    assert set(combo_iter) == combo_contacts

    duplicated = compute_contacts(topo, provider="molstar", only=[hb, hb, whb, filter_set])
    assert set(duplicated) == combo_contacts


def test_self_hydrophobic_like_count_and_invariants(luni: LahutaSystem, topo: Topology) -> None:
    sel = atoms(lambda a: a.type.has(AtomType.Hydrophobic))

    def within_45(i: int, j: int, d2: float) -> InteractionType:
        return InteractionType.Hydrophobic if d2 <= 4.5 * 4.5 else InteractionType.NoInteraction

    cs = find_contacts(topo, sel, tester=within_45, distance_max=4.5)
    assert cs.size() == 340

    keys = {_key(c) for c in cs}
    assert (Kind.Atom, 7, Kind.Atom, 9,  Category.Hydrophobic) in keys
    assert (Kind.Atom, 9, Kind.Atom, 10, Category.Hydrophobic) in keys

    hyp = pytest.importorskip("hypothesis",            reason="Hypothesis not installed")
    st  = pytest.importorskip("hypothesis.strategies", reason="Hypothesis not installed")

    @hyp.given(st.integers(min_value=0, max_value=max(0, cs.size() - 1)))
    def _check_item(i: int) -> None:
        c = cs[i]
        assert c.lhs.kind == Kind.Atom and c.rhs.kind == Kind.Atom
        assert c.type.category == Category.Hydrophobic
        assert c.distance_sq <= 4.5 * 4.5 + 1e-6

        # Both atoms should be hydrophobic by classification
        a_l = topo.get_atom(c.lhs.index)
        a_r = topo.get_atom(c.rhs.index)
        assert a_l.type.has(AtomType.Hydrophobic)
        assert a_r.type.has(AtomType.Hydrophobic)

    _check_item()


def test_cross_cation_pi_like_count_and_shape(luni: LahutaSystem, topo: Topology) -> None:
    sel_groups = groups(lambda g: g.a_type == AtomType.PositiveCharge or g.a_type == AtomType.NegativeCharge)
    sel_rings  = rings(lambda r: r.aromatic)

    def cation_pi(i: int, j: int, d2: float):
        return InteractionType.CationPi if d2 <= 6.0 * 6.0 else InteractionType.NoInteraction

    cs = find_contacts(topo, sel_groups, sel_rings, tester=cation_pi, distance_max=6.0)
    assert cs.size() == 7

    # Exemplars by indices (canonical Ring -> Group ordering)
    pair_keys = {_key(c) for c in cs if c.type.category == Category.CationPi}
    assert (Kind.Ring, 1, Kind.Group, 2,  Category.CationPi) in pair_keys
    assert (Kind.Ring, 5, Kind.Group, 7,  Category.CationPi) in pair_keys
    assert (Kind.Ring, 9, Kind.Group, 10, Category.CationPi) in pair_keys

    # Shape invariants (orientation may be canonicalized: Ring comes before Group)
    for c in cs:
        kinds = {c.lhs.kind, c.rhs.kind}
        assert kinds == {Kind.Ring, Kind.Group}
        if c.type.category != Category.Unclassified:
            assert c.type.category == Category.CationPi
            assert c.distance_sq >= 0.0


def test_self_pi_stacking_like_count(luni: LahutaSystem, topo: Topology) -> None:
    sel_rings = rings(lambda r: r.aromatic)

    def pi_stack(i: int, j: int, d2: float):
        return InteractionType.PiStacking if d2 <= 6.0 * 6.0 else InteractionType.NoInteraction

    cs = find_contacts(topo, sel_rings, tester=pi_stack, distance_max=6.0)
    assert cs.size() == 1

    c = cs[0]
    assert c.lhs.kind == Kind.Ring and c.rhs.kind == Kind.Ring
    assert c.type.category == Category.PiStacking
    assert c.distance_sq >= 0.0
