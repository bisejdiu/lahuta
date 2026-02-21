# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         def __init__(self): self._p = iter(["besian", "sejdiu", "@gmail.com"])
#         def __iter__(self): return self
#         def __next__(self): return next(self._p)
#     print("".join(Email()))
#
"""Using entities (EntityID/Kind) and records (AtomRec, RingRec, GroupRec) with RDKit atoms."""

from __future__ import annotations

from pathlib import Path

from lahuta import (
    AtomRec,
    AtomType,
    EntityID,
    GroupRec,
    Kind,
    LahutaSystem,
    RingRec,
    logging,
)
from lahuta.rdkit import Atom

DATA = Path(__file__).resolve().parents[3] / "data" / "ubi.cif"


# fmt: off
def main() -> None:
    sys = LahutaSystem(str(DATA))
    if not sys.build_topology():
        logging.error("Failed to build topology")
        return
    topo = sys.get_topology()

    # AtomRec and EntityID for an atom
    atoms: list[AtomRec] = topo.atom_records  # NOTE: type shown for clarity (it's correctly inferred)

    if not atoms:
        logging.warn("No AtomRec entries exposed")
        return

    a0 = atoms[0]
    atom_idx = int(a0.idx)
    logging.info(
        "AtomRec flags: donor=%s acceptor=%s hydrophobic=%s",
        a0.is_donor,
        a0.is_acceptor,
        a0.is_hydrophobic,
    )
    logging.info("Hydrophobic via contains: %s", AtomType.Hydrophobic in a0.type)
    logging.info("Hydrophobic atom count: %d", len(topo.atoms_with_type(AtomType.Hydrophobic)))

    eid_atom = EntityID.make(Kind.Atom, atom_idx)
    logging.info(f"Atom EntityID: {eid_atom} (repr={eid_atom!r}) kind={eid_atom.kind} index={eid_atom.index}")

    mol = topo.molecule()
    rd_atom = mol.getAtomWithIdx(atom_idx)
    logging.info(f"RDKit Atom[{atom_idx}]: Z={rd_atom.getAtomicNum()} symbol={rd_atom.getSymbol()} degree={rd_atom.getDegree()}")

    rings: list[RingRec] = topo.rings
    if rings:
        r0 = rings[0]
        atoms_in_ring: list[Atom] = r0.atoms

        logging.info(
            f"Ring[0]: size={r0.size} aromatic={r0.aromatic} "
            f"center=({r0.center.x:.3f},{r0.center.y:.3f},{r0.center.z:.3f}) "
            f"normal=({r0.normal.x:.3f},{r0.normal.y:.3f},{r0.normal.z:.3f})"
        )

        if atoms_in_ring and len(atoms_in_ring) >= 2:
            a0, a1 = atoms_in_ring[0], atoms_in_ring[1]
            logging.info(f"  first  ring atom: idx={a0.getIdx()} symbol={a0.getSymbol()} aromatic={a0.getIsAromatic()}")
            logging.info(f"  second ring atom: idx={a1.getIdx()} symbol={a1.getSymbol()} aromatic={a1.getIsAromatic()}")

        eid_ring = EntityID.make(Kind.Ring, 0)
        logging.info(f"Ring EntityID: {eid_ring}")

    groups: list[GroupRec] = topo.groups
    if groups:
        g0 = groups[0]
        atoms_in_group: list[Atom] = g0.atoms
        logging.info(f"Group[0]: n_atoms={len(atoms_in_group)} center=({g0.center.x:.3f},{g0.center.y:.3f},{g0.center.z:.3f})")
        if atoms_in_group: # groups can have just 1 atom
            a = atoms_in_group[0]
            logging.info(f"  first group atom: idx={a.getIdx()} symbol={a.getSymbol()}")

        eid_group = EntityID.make(Kind.Group, 0)
        logging.info(f"Group EntityID: {eid_group}")


if __name__ == "__main__":
    main()
