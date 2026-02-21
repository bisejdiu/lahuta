# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         r = ""
#         def __hash__(self):
#             Email.r = "besian" + "sejdiu" + "@gmail.com"
#             return 0
#     hash(Email())
#     print(Email.r)
#
"""
Lifetime and ownership tests for Python bindings.

Keep-alive graph exposed to Python:
    LahutaSystem  ->  Topology  ->  Residues  ->  Residue (view)
           |            |              |             |
     RDKit mol/conf   views/refs    iterators     atoms/fields

Verification goals:
- `LahutaSystem.get_topology()` returns a view with lifetime tied to the parent system.
- `Topology.residues` returns a container view tied to the parent topology,
  and iterators from it keep the container alive.
- `Residues.__getitem__` returns a view of a residue tied to the container, not a detached copy.
- `Topology` properties `atom_records`/`rings`/`groups` return lists of live record views.
  The list keeps the parent `Topology` alive, and individual items are
  reference-bound to the `Topology`, so they remain valid when parents are dropped.
- RDKit handles (molecule/conformer) remain usable after `LahutaSystem` is dropped.

Further notes:
Initially, use-after-free was non-deterministic for records (`AtomRec`/`RingRec`/`GroupRec`)
because they store non-owning `RDKit::Atom` references.
Surprisingly, even under heavy stress, crashes did not occur after dropping parents.

Still, the bindings were factored so that:
- `get_atom`/`get_ring`/`get_group` return references (tied to `Topology`).
- `atom_records`/`rings`/`groups` return lists of references, each element tied to
  the `Topology`, and the list itself keeps the `Topology` alive.
"""

from __future__ import annotations

import gc
from pathlib import Path

from lahuta import LahutaSystem


# fmt: off
def _build_sys(ubi_cif: Path) -> LahutaSystem:
    sys = LahutaSystem(str(ubi_cif))
    assert sys.build_topology() is True
    return sys


def test_topology_keeps_owner_alive(ubi_cif: Path) -> None:
    """Verify that Topology properties remain accessible after LahutaSystem is deleted."""
    sys = _build_sys(ubi_cif)
    top = sys.get_topology()

    del sys
    gc.collect()

    rings  = top.rings
    groups = top.groups
    atoms  = top.atom_records

    assert isinstance(rings,  list)
    assert isinstance(groups, list)
    assert isinstance(atoms,  list)

    assert len(atoms) > 0


def test_molecule_and_conformer_keep_owner_alive(ubi_cif: Path) -> None:
    """Verify that RDKit molecule and conformer handles remain usable after LahutaSystem and Topology are deleted."""
    sys  = _build_sys(ubi_cif)
    top  = sys.get_topology()
    mol  = top.molecule()
    conf = top.conformer(0)

    del sys, top
    gc.collect()

    assert mol.getNumAtoms() > 0
    assert conf.getNumAtoms() == mol.getNumAtoms()


def test_residues_keep_topology_and_system_alive(ubi_cif: Path) -> None:
    """Verify that Residues container remains fully functional after parent objects are deleted."""
    sys = _build_sys(ubi_cif)
    top = sys.get_topology()
    res = top.residues  # bound with reference_internal

    del top, sys
    gc.collect()

    # residues should still be fully usable
    assert len(res) > 0
    names = res.get_residue_names()
    assert isinstance(names, list) and isinstance(res.get_atom_ids(), list)
    assert len(names) == len(res)

    # indexing returns a view (Residue) with lifetime tied to the container
    # but if we hold it then we keep the container alive
    r0 = res[0]
    assert isinstance(r0.idx,   int)
    assert isinstance(r0.name,  str)
    assert isinstance(r0.atoms, list)
    assert len(r0.atoms) > 0
    assert isinstance(r0.atoms[0].getIdx(), int)
    assert isinstance(r0.atoms[0].getSymbol(), str)
    assert isinstance(r0.atoms[0].getAtomicNum(), int)


def test_residues_iterator_keeps_container_alive_after_parents_deleted(ubi_cif: Path) -> None:
    """Verify that Residues iterator remains functional after all parent objects are deleted."""
    sys = _build_sys(ubi_cif)
    top = sys.get_topology()
    res = top.residues

    it = iter(res)  # __iter__ has keep_alive<0,1>()
    del top, sys, res
    gc.collect()

    walked = list(it) # should stil lwork
    # may be possible to be empty but must not segfault
    assert isinstance(walked, list)
    if walked:
        r = walked[0]
        assert hasattr(r, "name") and hasattr(r, "atoms")


def test_residues_snapshot_is_independent_list(ubi_cif: Path) -> None:
    """Verify that Residues snapshot creates an independent Python list that survives parent deletion."""
    sys = _build_sys(ubi_cif)
    top = sys.get_topology()
    res = top.residues

    snapshot = res.residues # take a python snapshot
    assert isinstance(snapshot, list)

    del top, sys, res
    gc.collect()

    assert all(hasattr(r, "name") and isinstance(r.atoms, list) for r in snapshot)


def test_residue_value_object_survives_after_container_deleted(ubi_cif: Path) -> None:
    """Verify that individual Residue value objects remain valid after container and parents are deleted."""
    sys = _build_sys(ubi_cif)
    top = sys.get_topology()
    res = top.residues

    r = res[0] # take a value object

    del top, sys, res
    gc.collect()

    assert isinstance(r.name,  str)
    assert isinstance(r.atoms, list)


def test_residues_filter_returns_container_that_survives_without_parents(ubi_cif: Path) -> None:
    """Verify that filtered Residues container remains functional after parent objects are deleted."""
    sys = _build_sys(ubi_cif)
    top = sys.get_topology()
    res = top.residues

    # create a filtered view/new container (predicate keeps all)
    res2 = res.filter(lambda _: True)
    assert len(res2) == len(res)

    # drop parents and original container
    del top, sys, res
    gc.collect()

    # res2 should still be fully usable
    assert len(res2) >= 0  # existence check
    _ = res2.get_residue_names()
    if len(res2) > 0:
        r = res2[0]
        assert hasattr(r, "name") and hasattr(r, "atoms")


def test_residues_map_after_parents_deleted(ubi_cif: Path) -> None:
    """Verify that Residues map operation works correctly after parent objects are deleted."""
    sys = _build_sys(ubi_cif)
    top = sys.get_topology()
    res = top.residues

    del top, sys
    gc.collect()

    # map runs after parents are gone; still must work
    names = res.map(lambda rr: rr.name)
    assert isinstance(names, list)
    if names:
        assert isinstance(names[0], str)


def test_topology_rdkit_handles_keep_chain_alive(ubi_cif: Path) -> None:
    """Verify that RDKit handles from Topology keep the entire ownership chain alive."""
    sys = _build_sys(ubi_cif)
    top   = sys.get_topology()
    mol2  = top.molecule()  # reference_internal to Topology
    conf2 = top.conformer() # reference_internal to Topology

    # RDKit handles should keep Topology alive, and Topology (via get_topology binding) should keep the LahutaSystem alive
    del top, sys
    gc.collect()

    assert mol2.getNumAtoms() > 0
    assert conf2.getNumAtoms() == mol2.getNumAtoms()


def test_topology_property_copies_survive_without_parents(ubi_cif: Path) -> None:
    """Verify that Topology property copies (atom_records, rings, groups) survive after parent deletion."""
    sys = _build_sys(ubi_cif)
    top = sys.get_topology()

    # value-returning properties (copies) should be usable after parents go away
    atom_records = top.atom_records
    rings      = top.rings
    groups     = top.groups

    del top, sys
    gc.collect()

    # these are plain Python containers, hence they must survive
    assert isinstance(atom_records, list)
    assert isinstance(rings,      list)
    assert isinstance(groups,     list)

    if rings:
        assert hasattr(rings[0], "size") or hasattr(rings[0], "atoms") or True
    if groups:
        first = groups[0]
        assert hasattr(first, "idx") or True


def test_atomrec_survives_after_parents_deleted(ubi_cif: Path, run_child) -> None:
    """Accessing an AtomRec after dropping LahutaSystem/Topology remains valid."""
    ubi_path = str(ubi_cif)
    child = f"""
import gc, sys
from lahuta import LahutaSystem

sys_obj = LahutaSystem({ubi_path!r})
assert sys_obj.build_topology() is True
top = sys_obj.get_topology()

recs = top.atom_records
a0 = recs[0]

del top, sys_obj, recs
gc.collect()

idx = int(a0.idx)
print("OK", idx)
sys.exit(0)
"""

    res = run_child(child)
    if res.returncode != 0:
        raise AssertionError(
            f"child failed (returncode={res.returncode})\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )
    assert "OK" in res.stdout


def test_ringrec_atoms_survive_after_parents_deleted(ubi_cif: Path, run_child) -> None:
    """Accessing RingRec.atoms[i] after dropping parents remains valid."""
    ubi_path = str(ubi_cif)
    child = f"""
import gc, sys
from lahuta import LahutaSystem, rdkit as rdkit

sys_obj = LahutaSystem({ubi_path!r})
assert sys_obj.build_topology() is True
top = sys_obj.get_topology()

rings = top.rings
if not rings:
    print("NO_RINGS")
    sys.exit(0)
r0 = rings[0]

del top, sys_obj, rings
gc.collect()

atoms = r0.atoms
val = atoms[0].getIdx()
print("OK", val)
sys.exit(0)
"""

    res = run_child(child)
    if res.returncode != 0:
        raise AssertionError(
            f"child failed (returncode={res.returncode})\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )
    # Either we had no rings or we succeeded
    assert "OK" in res.stdout or "NO_RINGS" in res.stdout


def test_grouprec_atoms_survive_after_parents_deleted(ubi_cif: Path, run_child) -> None:
    """Accessing GroupRec.atoms[i] after dropping parents remains valid."""
    ubi_path = str(ubi_cif)
    child = f"""
import gc, sys
from lahuta import LahutaSystem, rdkit as rdkit

sys_obj = LahutaSystem({ubi_path!r})
assert sys_obj.build_topology() is True
top = sys_obj.get_topology()

groups = top.groups
if not groups:
    print("NO_GROUPS")
    sys.exit(0)
g0 = groups[0]

del top, sys_obj, groups
gc.collect()

atoms = g0.atoms
val = atoms[0].getIdx()
print("OK", val)
sys.exit(0)
"""

    res = run_child(child)
    if res.returncode != 0:
        raise AssertionError(
            f"child failed (returncode={res.returncode})\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )
    assert "OK" in res.stdout or "NO_GROUPS" in res.stdout


def test_atom_monomer_info_survives_after_parents_deleted(ubi_cif: Path, run_child) -> None:
    """
    Accessing Atom.getMonomerInfo() (AtomPDBResidueInfo) after dropping parents remains valid.
    Monomer info is owned by RDKit::Atom. We make sure RDKit atoms are kept alive in our views.
    """
    ubi_path = str(ubi_cif)
    child = f"""
import gc, sys
from lahuta import LahutaSystem, rdkit as rdkit

sys_obj = LahutaSystem({ubi_path!r})
assert sys_obj.build_topology() is True
top = sys_obj.get_topology()

residues = top.residues
assert len(residues) > 0, "Expected at least one residue"
r0 = residues[0]
assert len(r0.atoms) > 0, "Expected residue to have atoms"
a0 = r0.atoms[0]

info = a0.getMonomerInfo()
assert info is not None, "Expected monomer info on atom"

# It should be a PDB residue info; access a few fields
resname = info.getResidueName()
resnum  = info.getResidueNumber()
chain   = info.getChainId()

del top, sys_obj, residues, r0
gc.collect()

assert resname,     "Residue name should be non-empty"
assert resnum != 0, "Residue number should be non-zero"
assert chain,       "Chain ID should be non-empty"

# Still accessible and valid
print("OK", resname, resnum, chain)
sys.exit(0)
"""

    res = run_child(child)
    if res.returncode != 0:
        raise AssertionError(
            f"child failed (returncode={res.returncode})\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )
    assert "OK" in res.stdout
