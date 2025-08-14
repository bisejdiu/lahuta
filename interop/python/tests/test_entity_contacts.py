"""Tests for entity-based contact detection."""

from __future__ import annotations

import pytest

import lahuta as lxx
from lahuta.entities import atoms, compute_contacts, find_contacts, groups, rings


# fmt: off
@pytest.fixture(scope="session")
def topo(luni: lxx.LahutaSystem) -> lxx.Topology:
    assert luni.build_topology() is True
    return luni.get_topology()


def _key(c: lxx.Contact) -> tuple[lxx.Kind, int, lxx.Kind, int, lxx.Category]:
    d = c.to_dict()
    return d["lhs_kind"], int(d["lhs_index"]), d["rhs_kind"], int(d["rhs_index"]), d["category"]


def test_system_sanity(luni: lxx.LahutaSystem, topo: lxx.Topology) -> None:
    assert isinstance(luni.n_atoms, int) and luni.n_atoms == 671
    assert len(topo.atom_types) == 671
    assert len(topo.rings)  >= 1
    assert len(topo.groups) >= 1


def test_provider_engine_molstar_counts_and_samples(luni: lxx.LahutaSystem, topo: lxx.Topology) -> None:
    ms_all = compute_contacts(topo, provider="molstar")
    hbonds = compute_contacts(topo, provider="molstar", only=lxx.InteractionType.HydrogenBond)

    assert ms_all.size() == 153
    assert hbonds.size() == 99

    # Spot-check exemplars via indices and category using structured data
    keys_all = {_key(c) for c in ms_all}
    keys_hb  = {_key(c) for c in hbonds}

    def atom_pair(i: int, j: int, cat: lxx.Category) -> tuple[lxx.Kind, int, lxx.Kind, int, lxx.Category]:
        return (lxx.Kind.Atom, min(i, j), lxx.Kind.Atom, max(i, j), cat)

    assert atom_pair(4,  76, lxx.Category.HydrogenBond) in keys_all
    assert atom_pair(8,  72, lxx.Category.HydrogenBond) in keys_all
    assert atom_pair(28, 52, lxx.Category.WeakHydrogenBond) in keys_all
    assert atom_pair(29, 37, lxx.Category.HydrogenBond) in keys_all

    assert atom_pair(4, 76,  lxx.Category.HydrogenBond) in keys_hb
    assert atom_pair(8, 72,  lxx.Category.HydrogenBond) in keys_hb
    assert atom_pair(29, 37, lxx.Category.HydrogenBond) in keys_hb


def test_self_hydrophobic_like_count_and_invariants(luni: lxx.LahutaSystem, topo: lxx.Topology) -> None:
    sel = atoms(lambda a: a.type.has(lxx.AtomType.Hydrophobic))

    def within_45(i: int, j: int, d2: float) -> lxx.InteractionType:
        return lxx.InteractionType.Hydrophobic if d2 <= 4.5 * 4.5 else lxx.InteractionType.None_

    cs = find_contacts(topo, sel, tester=within_45, distance_max=4.5)
    assert cs.size() == 340

    keys = {_key(c) for c in cs}
    assert (lxx.Kind.Atom, 7, lxx.Kind.Atom, 9,  lxx.Category.Hydrophobic) in keys
    assert (lxx.Kind.Atom, 9, lxx.Kind.Atom, 10, lxx.Category.Hydrophobic) in keys

    hyp = pytest.importorskip("hypothesis",            reason="Hypothesis not installed")
    st  = pytest.importorskip("hypothesis.strategies", reason="Hypothesis not installed")

    @hyp.given(st.integers(min_value=0, max_value=max(0, cs.size() - 1)))
    def _check_item(i: int) -> None:
        c = cs[i]
        assert c.lhs.kind == lxx.Kind.Atom and c.rhs.kind == lxx.Kind.Atom
        assert c.type.category == lxx.Category.Hydrophobic
        assert c.distance_sq <= 4.5 * 4.5 + 1e-6

        # Both atoms should be hydrophobic by classification
        a_l = topo.get_atom(c.lhs.index)
        a_r = topo.get_atom(c.rhs.index)
        assert a_l.type.has(lxx.AtomType.Hydrophobic)
        assert a_r.type.has(lxx.AtomType.Hydrophobic)

    _check_item()


def test_cross_cation_pi_like_count_and_shape(luni: lxx.LahutaSystem, topo: lxx.Topology) -> None:
    sel_groups = groups(lambda g: g.a_type == lxx.AtomType.PositiveCharge or g.a_type == lxx.AtomType.NegativeCharge)
    sel_rings  = rings(lambda r: r.aromatic)

    def cation_pi(i: int, j: int, d2: float):
        return lxx.InteractionType.CationPi if d2 <= 6.0 * 6.0 else lxx.InteractionType.None_

    cs = find_contacts(topo, sel_groups, sel_rings, tester=cation_pi, distance_max=6.0)
    assert cs.size() == 7

    # Exemplars by indices (canonical Ring -> Group ordering)
    pair_keys = {_key(c) for c in cs if c.type.category == lxx.Category.CationPi}
    assert (lxx.Kind.Ring, 1, lxx.Kind.Group, 2,  lxx.Category.CationPi) in pair_keys
    assert (lxx.Kind.Ring, 5, lxx.Kind.Group, 7,  lxx.Category.CationPi) in pair_keys
    assert (lxx.Kind.Ring, 9, lxx.Kind.Group, 10, lxx.Category.CationPi) in pair_keys

    # Shape invariants (orientation may be canonicalized: Ring comes before Group)
    for c in cs:
        kinds = {c.lhs.kind, c.rhs.kind}
        assert kinds == {lxx.Kind.Ring, lxx.Kind.Group}
        if c.type.category != lxx.Category.None_:
            assert c.type.category == lxx.Category.CationPi
            assert c.distance_sq >= 0.0


def test_self_pi_stacking_like_count(luni: lxx.LahutaSystem, topo: lxx.Topology) -> None:
    sel_rings = rings(lambda r: r.aromatic)

    def pi_stack(i: int, j: int, d2: float):
        return lxx.InteractionType.PiStacking if d2 <= 6.0 * 6.0 else lxx.InteractionType.None_

    cs = find_contacts(topo, sel_rings, tester=pi_stack, distance_max=6.0)
    assert cs.size() == 1

    c = cs[0]
    assert c.lhs.kind == lxx.Kind.Ring and c.rhs.kind == lxx.Kind.Ring
    assert c.type.category == lxx.Category.PiStacking
    assert c.distance_sq >= 0.0
