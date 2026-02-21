# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     s = ""
#     s += "besian"
#     s += "sejdiu"
#     s += "@gmail.com"
#     print(s)
#
"""
Resolving EntityID to concrete records from Python.
Prefer Contact.describe(topology) / Contact.to_dict(topology) for quick inspection.
See EntityResolver for more.
"""

from __future__ import annotations

from pathlib import Path

from lahuta import EntityResolver, LahutaSystem, MolStarContactsEngine, Topology
from lahuta.lib.lahuta import ContactSet

DATA = Path(__file__).resolve().parents[3] / "core" / "data" / "ubi.cif"


# fmt: off
def build_topology_from_file() -> Topology:
    luni = LahutaSystem(str(DATA))
    assert luni.build_topology(), "Failed to build topology"
    return luni.get_topology()


def demo_contact_descriptions(top: Topology, cset: ContactSet) -> None:
    """The recommended helpers for contacts."""
    print("Contact.describe (first 3 contacts)")
    for i in range(min(3, cset.size())):
        print("  ", cset[i].describe(top))

    if cset.size():
        print("Contact.to_dict (first contact)")
        print("  ", cset[0].to_dict(top))


def demo_entity_resolver_minimal(top: Topology, cset: ContactSet) -> None:
    """Resolve EntityID to concrete records when needed."""
    if cset.empty():
        return
    resolver = EntityResolver(top)
    lhs_rec, rhs_rec = resolver.resolve_contact(cset[0])
    print("EntityResolver.resolve_contact (first contact)")
    print("  ", lhs_rec, "-&-", rhs_rec)


def main() -> None:
    top = build_topology_from_file()

    engine = MolStarContactsEngine()
    cset = engine.compute(top)
    print(f"Contacts: {cset.size()}")

    demo_contact_descriptions(top, cset)
    demo_entity_resolver_minimal(top, cset)


if __name__ == "__main__":
    main()
