"""Topology bindings: flags, options, residues container, typing assignment, and basic accessors."""

from __future__ import annotations

from pathlib import Path

import pytest

import lahuta as lxx


# Isolated, per-test system with topology already built.
# Avoids test ordering dependence when running with pytest-xdist.
@pytest.fixture(scope="function")
def luni_built(ubi_cif: Path) -> lxx.LahutaSystem:
    sys = lxx.LahutaSystem(str(ubi_cif))
    assert sys.build_topology() is True
    return sys


# fmt: off
def test_topology_build_with_options_and_flags(luni: lxx.LahutaSystem) -> None:
    # Configure explicit options
    opts = lxx.TopologyBuildingOptions()
    opts.cutoff = 4.5
    opts.compute_nonstandard_bonds = True
    opts.atom_typing_method = lxx.AtomTypingMethod.Molstar

    # Start with no stages, then selectively enable
    luni.enable_only(lxx.TopologyComputers.None_)
    luni.enable_only(lxx.TopologyComputers.Standard)

    assert luni.build_topology(opts) is True
    assert luni.has_topology_built() is True

    if not luni.is_computation_enabled(lxx.TopologyComputers.Rings):
        luni.enable_computation(lxx.TopologyComputers.Rings, True)
    assert luni.execute_computation(lxx.TopologyComputers.Rings) is True


def test_residues_container_and_helpers(luni_built: lxx.LahutaSystem) -> None:
    topo = luni_built.get_topology()
    residues = topo.residues

    n = len(residues)
    assert n > 0

    # Names via API should match map over Residue.name
    names_api = residues.get_residue_names()
    names_map = residues.map(lambda r: r.name)
    assert isinstance(names_api, list) and isinstance(names_map, list)
    assert len(names_api) == n == len(names_map)

    # Filter returns non-empty for any known name
    sub = residues.filter(lambda r: r.name == names_api[0])
    assert len(sub) >= 1

    # Atom IDs is a flat list of indices
    atom_ids = residues.get_atom_ids()
    assert isinstance(atom_ids, list)
    assert all(isinstance(i, int) and i >= 0 for i in atom_ids)


def test_atom_typing_and_records(luni_built: lxx.LahutaSystem) -> None:
    topo = luni_built.get_topology()

    # Assign types using both backends, lists must have size N_atoms
    topo.set_atom_typing_method(lxx.AtomTypingMethod.Molstar)
    topo.assign_typing(lxx.AtomTypingMethod.Molstar)
    types_molstar = topo.atom_types

    topo.set_atom_typing_method(lxx.AtomTypingMethod.Arpeggio)
    topo.assign_typing(lxx.AtomTypingMethod.Arpeggio)
    types_arpeggio = topo.atom_types

    assert isinstance(types_molstar, list) and isinstance(types_arpeggio, list)
    assert len(types_molstar) == len(types_arpeggio) == luni_built.n_atoms

    # Filter-by-fn should preserve count when predicate is True for all
    all_recs = topo.atom_types_filter_by_fn(lambda _: True)
    assert isinstance(all_recs, list) and len(all_recs) == luni_built.n_atoms

    # Access a single atom record and check idx consistency
    rec0 = topo.get_atom(0)
    assert hasattr(rec0, "idx")
    assert isinstance(rec0.idx(), int) and rec0.idx() == 0


def test_molecule_and_conformer_handles(luni_built: lxx.LahutaSystem) -> None:
    topo = luni_built.get_topology()
    mol  = topo.molecule()
    conf = topo.conformer()
    xyz  = luni_built.props.positions

    assert mol.getNumAtoms()  == luni_built.n_atoms
    assert conf.getNumAtoms() == luni_built.n_atoms
    assert xyz.shape == (luni_built.n_atoms, 3)
