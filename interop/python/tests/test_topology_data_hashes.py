"""Strict hashing tests for topology data: atoms, rings, groups."""

from __future__ import annotations

import hashlib
import struct
from pathlib import Path

import pytest

import lahuta as lxx


# fmt: off
def _sha() -> "hashlib._Hash":
    return hashlib.sha256()

def hash_atoms(atom_recs) -> str:
    h = _sha()
    # header: tag + count
    h.update(b"LHX2|atoms|")
    h.update(struct.pack("<Q", len(atom_recs)))
    for rec in atom_recs:
        idx = int(rec.idx())
        at = int(rec.type)
        h.update(struct.pack("<II", idx, at))
    return h.hexdigest()

def _pack_point3d(h, pt) -> None:
    h.update(struct.pack("<ddd", float(pt.x), float(pt.y), float(pt.z)))

def hash_rings(ring_recs) -> str:
    h = _sha()
    h.update(b"LHX2|rings|")
    h.update(struct.pack("<Q", len(ring_recs)))
    for r in ring_recs:
        atoms = r.atoms  # list of RDKit Atom
        n = len(atoms)
        h.update(struct.pack("<I", n))
        for a in atoms:
            # RDKit Python Atom exposes getIdx()
            h.update(struct.pack("<I", int(a.getIdx())))
        h.update(struct.pack("<?", bool(r.aromatic)))
        _pack_point3d(h, r.center)
        _pack_point3d(h, r.normal)
    return h.hexdigest()

def hash_groups(group_recs) -> str:
    h = _sha()
    h.update(b"LHX2|groups|")
    h.update(struct.pack("<Q", len(group_recs)))
    for g in group_recs:
        fg = int(g.type)
        at = int(g.a_type)
        h.update(struct.pack("<II", fg, at))
        atoms = g.atoms  # list of RDKit Atom
        n = len(atoms)
        h.update(struct.pack("<I", n))
        for a in atoms:
            h.update(struct.pack("<I", int(a.getIdx())))
        _pack_point3d(h, g.center)
    return h.hexdigest()


EXPECTED = {
    "atoms": {
        "count": 671,
        "sha256": "d20618dea32021ec300dc1985a38388fc4db677c10fc64504f2151137ba031fb",
    },
    "rings": {
        "count": 10,
        "sha256": "69ef8c86f68085e5ece89e5a0f865a46058d779e5105673e323be8be4470be3a",
    },
    "groups": {
        "count": 35,
        "sha256": "ef8a0dfc5b456566fb797f51a58592b9aba16bba199711d5b481d27fb6679a6a",
    },
}


@pytest.fixture(scope="session")
def luni(ubi_cif: Path) -> lxx.LahutaSystem:
    return lxx.LahutaSystem(str(ubi_cif))


def _ensure_topology(luni: lxx.LahutaSystem) -> None:
    # Build with a stable stage selection
    opts = lxx.TopologyBuildingOptions()
    opts.cutoff = 4.5
    assert luni.build_topology(opts) is True


def test_topology_hashes_atoms_rings_groups(luni: lxx.LahutaSystem, capsys: pytest.CaptureFixture[str]) -> None:
    _ensure_topology(luni)
    topo = luni.get_topology()

    atoms  = topo.atom_types
    rings  = topo.rings
    groups = topo.groups

    got = {
        "atoms": {
            "count": len(atoms),
            "sha256": hash_atoms(atoms),
        },
        "rings": {
            "count": len(rings),
            "sha256": hash_rings(rings),
        },
        "groups": {
            "count": len(groups),
            "sha256": hash_groups(groups),
        },
    }

    print("Computed topology digests:")
    for k in ("atoms", "rings", "groups"):
        print(f"  {k}: count={got[k]['count']} sha256={got[k]['sha256']}")
    captured = capsys.readouterr()
    _ = captured.out

    for k in ("atoms", "rings", "groups"):
        exp = EXPECTED[k]
        if exp["sha256"].startswith("<fill"):
            pytest.xfail("Update EXPECTED digests in test_topology_data_hashes.py with printed values.")

    for k in ("atoms", "rings", "groups"):
        assert got[k]["count"]  == EXPECTED[k]["count"],  f"{k}: count mismatch"
        assert got[k]["sha256"] == EXPECTED[k]["sha256"], f"{k}: sha256 mismatch"
