# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     def f(**kw):
#         return kw["f"] + kw["l"] + "@" + kw["d"]
#     print(f(**{"f": "besian", "l": "sejdiu", "d": "gmail.com"}))
#
"""Topology bindings: flags, options, residues container, typing assignment, and basic accessors."""

from __future__ import annotations

from pathlib import Path

import pytest

from lahuta import AtomType, AtomTypingMethod, LahutaSystem, TopologyBuildingOptions, TopologyComputers


# Isolated, per-test system with topology already built.
# Avoids test ordering dependence when running with pytest-xdist.
@pytest.fixture(scope="function")
def luni_built(ubi_cif: Path) -> LahutaSystem:
    sys = LahutaSystem(str(ubi_cif))
    assert sys.build_topology() is True
    return sys


# fmt: off
def test_topology_build_with_options_and_flags(luni: LahutaSystem) -> None:
    # Configure explicit options
    opts = TopologyBuildingOptions()
    opts.cutoff = 4.5
    opts.compute_nonstandard_bonds = True
    opts.atom_typing_method = AtomTypingMethod.MolStar

    # Build standard topology first
    assert luni.build_topology(opts, include=TopologyComputers.Standard) is True
    assert luni.has_topology_built() is True

    # Ensure rings are computed
    assert luni.build_topology(opts, include=TopologyComputers.Rings) is True
    topo = luni.get_topology()
    assert topo.has_computed(TopologyComputers.Rings) is True


def test_build_topology_include_only_defaults(ubi_cif: Path) -> None:
    sys = LahutaSystem(str(ubi_cif))
    assert sys.build_topology(include=TopologyComputers.Neighbors) is True
    topo = sys.get_topology()
    assert topo.has_computed(TopologyComputers.Neighbors) is True
    assert topo.has_computed(TopologyComputers.Bonds) is False


def test_build_topology_include_list_defaults(ubi_cif: Path) -> None:
    sys = LahutaSystem(str(ubi_cif))
    assert sys.build_topology(include=[TopologyComputers.Residues]) is True
    topo = sys.get_topology()
    assert topo.has_computed(TopologyComputers.Residues) is True
    assert topo.has_computed(TopologyComputers.Bonds) is False


def test_system_residues_lazy(ubi_cif: Path) -> None:
    sys = LahutaSystem(str(ubi_cif))
    residues = sys.residues
    assert len(residues) > 0
    assert sys.has_topology_built() is False


def test_get_or_build_topology_defaults(ubi_cif: Path) -> None:
    sys = LahutaSystem(str(ubi_cif))
    topo = sys.get_or_build_topology()
    assert topo is not None
    assert sys.has_topology_built() is True
    assert topo.has_computed(TopologyComputers.Standard) is True


def test_get_or_build_topology_with_options(ubi_cif: Path) -> None:
    sys = LahutaSystem(str(ubi_cif))
    opts = TopologyBuildingOptions()
    opts.cutoff = 3.9
    topo = sys.get_or_build_topology(opts)
    assert topo is not None
    assert sys.has_topology_built() is True


def test_reset_topology_file_backed(ubi_cif: Path) -> None:
    sys = LahutaSystem(str(ubi_cif))
    assert sys.build_topology(include=TopologyComputers.Neighbors) is True

    fresh = sys.reset_topology()
    assert fresh is not sys
    assert fresh.n_atoms == sys.n_atoms

    assert fresh.build_topology(include=TopologyComputers.Neighbors) is True
    assert sys.build_topology(include=TopologyComputers.Bonds) is True


def test_reset_topology_non_file_backed_raises(ubi_cif: Path) -> None:
    sys = LahutaSystem(str(ubi_cif))
    filtered = sys.filter([0, 1, 2])
    with pytest.raises(RuntimeError):
        filtered.reset_topology()


def test_residues_container_and_helpers(luni_built: LahutaSystem) -> None:
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


def test_atom_typing_and_records(luni_built: LahutaSystem) -> None:
    topo = luni_built.get_topology()

    # Assign types using both backends, lists must have size N_atoms
    topo.assign_typing(AtomTypingMethod.MolStar)
    types_molstar = topo.atom_records

    topo.assign_typing(AtomTypingMethod.Arpeggio)
    types_arpeggio = topo.atom_records

    assert isinstance(types_molstar, list) and isinstance(types_arpeggio, list)
    assert len(types_molstar) == len(types_arpeggio) == luni_built.n_atoms

    # Filter-by-fn should preserve count when predicate is True for all
    all_recs = topo.atom_types_filter_by_fn(lambda _: True)
    assert isinstance(all_recs, list) and len(all_recs) == luni_built.n_atoms

    # atoms_with_type should match manual filtering
    hydrophobic = topo.atoms_with_type(AtomType.Hydrophobic)
    manual = [r for r in topo.atom_records if AtomType.Hydrophobic in r.type]
    assert isinstance(hydrophobic, list)
    assert len(hydrophobic) == len(manual)

    # Access a single atom record and check idx consistency
    rec0 = topo.get_atom(0)
    assert hasattr(rec0, "idx")
    assert isinstance(rec0.idx, int) and rec0.idx == 0


def test_molecule_and_conformer_handles(luni_built: LahutaSystem) -> None:
    topo = luni_built.get_topology()
    mol  = topo.molecule()
    conf = topo.conformer(0)
    xyz  = luni_built.props.positions

    assert mol.getNumAtoms()  == luni_built.n_atoms
    assert conf.getNumAtoms() == luni_built.n_atoms
    assert xyz.shape == (luni_built.n_atoms, 3)
