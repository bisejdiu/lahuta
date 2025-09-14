"""
This example demonstrates the high-level contacts API for a few common cases.
The goal is to showcase the API surface. It does not reproduce the exact provider implementations.
"""

from __future__ import annotations

from pathlib import Path

import lahuta as lxx
from lahuta.entities import atoms, compute_contacts, find_contacts, groups, rings

DATA = Path(__file__).resolve().parents[3] / "data" / "ubi.cif"


# fmt: off
def build_topology() -> tuple[lxx.LahutaSystem, lxx.Topology]:
    sys = lxx.LahutaSystem(str(DATA))

    if not sys.build_topology():
        raise RuntimeError("Failed to build topology for contact examples")
    return sys, sys.get_topology()


def example_provider_engine(top: lxx.Topology) -> None:
    all_cs = compute_contacts(top, provider="molstar")
    print(f"MolStar: all interactions = {all_cs.size()}")

    hb = compute_contacts(top, provider="molstar", only=lxx.InteractionType.HydrogenBond)
    print(f"MolStar: hydrogen bonds = {hb.size()}")


def example_self_hydrophobic_like(top: lxx.Topology) -> None:
    # self-search: one custom filter + tester (hydrophobic-like)
    sel = atoms(lambda atom_rec: atom_rec.type.has(lxx.AtomType.Hydrophobic))

    def within_45_tester(i: int, j: int, d2: float) -> lxx.InteractionType:
        return lxx.InteractionType.Hydrophobic if d2 <= 4.5 * 4.5 else lxx.InteractionType.None_

    cs = find_contacts(top, sel, tester=within_45_tester, distance_max=4.5)
    print(f"Hydrophobic-like atom-atom contacts: {cs.size()}")


def example_cross_cation_pi_like(top: lxx.Topology) -> None:
    # cross-search: two custom filters + tester (cation-pi-like)
    sel_groups = groups(lambda g: g.a_type == lxx.AtomType.PositiveCharge or g.a_type == lxx.AtomType.NegativeCharge)
    sel_rings = rings(lambda r: r.aromatic)

    def cation_pi_tester(i: int, j: int, d2: float) -> lxx.InteractionType:
        return lxx.InteractionType.CationPi if d2 <= 6.0 * 6.0 else lxx.InteractionType.None_

    cs = find_contacts(top, sel_groups, sel_rings, tester=cation_pi_tester, distance_max=6.0)
    print(f"Cation-pi-like group-ring contacts: {cs.size()}")


def example_self_pi_stacking_like(top: lxx.Topology) -> None:
    # self-search: one custom filter + tester (pi-stacking-like)
    sel_rings = rings(lambda r: r.aromatic)

    def pi_stack_tester(i: int, j: int, d2: float) -> lxx.InteractionType:
        return lxx.InteractionType.PiStacking if d2 <= 6.0 * 6.0 else lxx.InteractionType.None_

    cs = find_contacts(top, sel_rings, tester=pi_stack_tester, distance_max=6.0)
    print(f"Pi-stacking-like ring-ring contacts: {cs.size()}")


def main() -> None:
    _, top = build_topology()
    example_provider_engine(top)
    example_self_hydrophobic_like(top)
    example_cross_cation_pi_like(top)
    example_self_pi_stacking_like(top)


if __name__ == "__main__":
    main()
