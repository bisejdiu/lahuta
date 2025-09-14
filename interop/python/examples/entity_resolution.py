"""
Resolving EntityID to concrete records from Python.
Demonstrates both Topology's typed resolving helpers and the EntityResolver layer.
"""

from __future__ import annotations

from pathlib import Path

from lahuta import EntityResolver, Kind, LahutaSystem, MolStarContactsEngine, Topology
from lahuta.lib.lahuta import ContactSet

DATA = Path(__file__).resolve().parents[3] / "data" / "ubi.cif"


# fmt: off
def build_topology_from_file() -> Topology:
    luni = LahutaSystem(str(DATA))
    assert luni.build_topology(), "Failed to build topology"
    return luni.get_topology()


def demo_topology_resolve(top: Topology, cset: ContactSet) -> None:
    """Show typed resolution using Topology.resolve_* for a few contacts."""
    print("Topology.resolve_* (first 3 contacts)")

    # NOTE: Type inference works correctly.
    for i in range(min(3, cset.size())):
        c = cset[i]
        if c.lhs.kind == Kind.Atom:
            lhs = top.resolve_atom(c.lhs)
        elif c.lhs.kind == Kind.Ring:
            lhs = top.resolve_ring(c.lhs)
        else:
            lhs = top.resolve_group(c.lhs)

        if c.rhs.kind == Kind.Atom:
            rhs = top.resolve_atom(c.rhs)
        elif c.rhs.kind == Kind.Ring:
            rhs = top.resolve_ring(c.rhs)
        else:
            rhs = top.resolve_group(c.rhs)

        print("  ", lhs, "-&-", rhs)


def demo_resolver_single(top: Topology, cset: ContactSet) -> None:
    """Use EntityResolver.resolve_contact on the first few contacts."""
    print("EntityResolver.resolve_contact (first 3 contacts)")
    resolver = EntityResolver(top)
    for i in range(min(3, cset.size())):
        lhs_rec, rhs_rec = resolver.resolve_contact(cset[i])
        print("  ", lhs_rec, "-&-", rhs_rec)


def demo_resolver_all(top: Topology, cset) -> None:
    """Materialize all pairs with EntityResolver.resolve_all and preview the first few."""
    print("EntityResolver.resolve_all (first 3 pairs)")
    resolver = EntityResolver(top)
    pairs = resolver.resolve_all(cset)
    for a, b in pairs[:3]:
        print("  ", a, "-&-", b)


def main() -> None:
    top = build_topology_from_file()

    engine = MolStarContactsEngine()
    cset = engine.compute(top)
    print(f"Contacts: {cset.size()}")

    demo_topology_resolve(top, cset)
    demo_resolver_single(top, cset)
    demo_resolver_all(top, cset)


if __name__ == "__main__":
    main()
