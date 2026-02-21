# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     print(str.__new__(str, "besian") + str.__new__(str, "sejdiu") + str.__new__(str, "@gmail.com"))
#
"""
This example demonstrates the high-level contacts API for a few common cases.
The goal is to showcase the API surface. It does not reproduce the exact provider implementations.
"""

from __future__ import annotations

from pathlib import Path

from lahuta import AtomType, InteractionType, LahutaSystem, Topology
from lahuta.entities import atoms, compute_contacts, find_contacts, groups, rings

DATA = Path(__file__).resolve().parents[3] / "data" / "ubi.cif"


# fmt: off
def build_topology() -> tuple[LahutaSystem, Topology]:
    sys = LahutaSystem(str(DATA))

    if not sys.build_topology():
        raise RuntimeError("Failed to build topology for contact examples")
    return sys, sys.get_topology()


def example_provider_engine(top: Topology) -> None:
    all_cs = compute_contacts(top, provider="molstar")
    print(f"MolStar: all interactions = {all_cs.size()}")

    hb = compute_contacts(top, provider="molstar", only=InteractionType.HydrogenBond)
    print(f"MolStar: hydrogen bonds = {hb.size()}")


def example_self_hydrophobic_like(top: Topology) -> None:
    # self-search: one custom filter + tester (hydrophobic-like)
    hydrophobic_atom_selector = atoms(lambda atom_rec: atom_rec.type.has(AtomType.Hydrophobic))

    def within_45_tester(i: int, j: int, d2: float) -> InteractionType:
        return InteractionType.Hydrophobic if d2 <= 4.5 * 4.5 else InteractionType.NoInteraction

    cs = find_contacts(top, hydrophobic_atom_selector, tester=within_45_tester, distance_max=4.5)
    print(f"Hydrophobic-like atom-atom contacts: {cs.size()}")


def example_cross_cation_pi_like(top: Topology) -> None:
    # cross-search: two custom filters + tester (cation-pi-like)
    charged_group_selector = groups(lambda g: g.a_type == AtomType.PositiveCharge or g.a_type == AtomType.NegativeCharge)
    aromatic_ring_selector = rings(lambda r: r.aromatic)

    def cation_pi_tester(i: int, j: int, d2: float) -> InteractionType:
        return InteractionType.CationPi if d2 <= 6.0 * 6.0 else InteractionType.NoInteraction

    cs = find_contacts(top, charged_group_selector, aromatic_ring_selector, tester=cation_pi_tester, distance_max=6.0)
    print(f"Cation-pi-like group-ring contacts: {cs.size()}")


def example_self_pi_stacking_like(top: Topology) -> None:
    # self-search: one custom filter + tester (pi-stacking-like)
    aromatic_ring_selector = rings(lambda r: r.aromatic)

    def pi_stack_tester(i: int, j: int, d2: float) -> InteractionType:
        return InteractionType.PiStacking if d2 <= 6.0 * 6.0 else InteractionType.NoInteraction

    cs = find_contacts(top, aromatic_ring_selector, tester=pi_stack_tester, distance_max=6.0)
    print(f"Pi-stacking-like ring-ring contacts: {cs.size()}")


def main() -> None:
    _, top = build_topology()
    example_provider_engine(top)
    example_self_hydrophobic_like(top)
    example_cross_cation_pi_like(top)
    example_self_pi_stacking_like(top)


if __name__ == "__main__":
    main()
