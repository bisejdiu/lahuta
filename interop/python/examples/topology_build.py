# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     print(operator.concat(operator.concat("besian", "sejdiu"), "@gmail.com"))
#
"""Controllig topology building and accessing topology data."""

from pathlib import Path
from typing import Callable

from lahuta import AtomTypingMethod, LahutaSystem, Residue, TopologyComputers, logging


# fmt: off
def topology_build(path: str | Path) -> None:
    logging.set_global_verbosity(logging.LogLevel.INFO)
    sys = LahutaSystem(str(path))
    sys.enable_only(TopologyComputers.Bonds)
    sys.set_search_cutoff_for_bonds(5)
    sys.build_topology()

    if not sys.is_computation_enabled(TopologyComputers.Rings):
        sys.enable_computation(TopologyComputers.Rings, True)
    sys.execute_computation(TopologyComputers.Rings)

    top = sys.get_topology()
    logging.info(f"Number of rings: {len(top.rings)}")
    logging.info("topology built")


def topology_counts_and_rings(path: str | Path) -> LahutaSystem:
    sys = LahutaSystem(str(path))
    sys.build_topology()
    topo = sys.get_topology()

    n_res    = len(topo.residues)
    n_groups = len(topo.groups)

    res_filter: Callable[[Residue], None] = lambda r: logging.info(f"Residue {r.name} ({r.idx})") if r.name.startswith("G") else None # noqa: E731
    topo.residues.map(res_filter)

    logging.info(f"Counts: residues={n_res}, groups={n_groups}")
    return sys

def typing_compare(sys: LahutaSystem) -> None:
    top = sys.get_topology()

    top.assign_typing(AtomTypingMethod.MolStar)
    mol_types = [rec.type for rec in top.atom_types]

    top.assign_typing(AtomTypingMethod.Arpeggio)
    arp_types = [rec.type for rec in top.atom_types]

    diffs = sum(m != a for m, a in zip(mol_types, arp_types))
    logging.info(f"typing differences: {diffs} / {len(mol_types)}")


if __name__ == "__main__":
    DATA = Path(__file__).resolve().parents[3] / "data" / "ubi.cif"

    topology_build(DATA)
    sys = topology_counts_and_rings(DATA)
    typing_compare(sys)
