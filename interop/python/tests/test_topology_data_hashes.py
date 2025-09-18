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

def _pack_point3d(h, pt, prec: int = 6) -> None:
    # Quantize to improve determinism
    h.update(
        struct.pack(
            "<ddd",
            round(float(pt.x), prec),
            round(float(pt.y), prec),
            round(float(pt.z), prec),
        )
    )

def _canonical_cycle(indices: list[int]) -> tuple[int, ...]:
    """Make ring atom ordering deterministic across diff platforms & compilers."""
    n = len(indices)
    if n == 0:
        return ()
    m = min(indices)
    positions = [i for i, v in enumerate(indices) if v == m]
    best: tuple[int, ...] | None = None
    for pos in positions:
        rot = tuple(indices[pos:] + indices[:pos])
        rrot = tuple(reversed(rot))
        cand = rot if rot < rrot else rrot
        if best is None or cand < best:
            best = cand
    assert best is not None
    return best

def hash_rings(ring_recs) -> str:
    entries = []
    for r in ring_recs:
        idxs = [int(a.getIdx()) for a in r.atoms]
        canon = _canonical_cycle(idxs)
        entries.append(
            {
                "idxs":     canon,
                "aromatic": bool(r.aromatic),
                "center":   r.center,
                "normal":   r.normal,
            }
        )
    entries.sort(key=lambda e: e["idxs"])

    h = _sha()
    h.update(b"LHX2|rings|")
    h.update(struct.pack("<Q", len(entries)))
    for e in entries:
        idxs = e["idxs"]
        h.update(struct.pack("<I", len(idxs)))
        for idx in idxs:
            h.update(struct.pack("<I", idx))
        h.update(struct.pack("<?", e["aromatic"]))
        _pack_point3d(h, e["center"], prec=6)
        _pack_point3d(h, e["normal"], prec=6)
    return h.hexdigest()

def hash_groups(group_recs) -> str:
    entries = []
    for g in group_recs:
        fg = int(g.type)
        at = int(g.a_type)
        idxs = [int(a.getIdx()) for a in g.atoms]
        canon = _canonical_cycle(idxs)
        entries.append({"fg": fg, "at": at, "idxs": canon, "center": g.center})
    entries.sort(key=lambda e: (e["fg"], e["at"], e["idxs"]))

    h = _sha()
    h.update(b"LHX2|groups|")
    h.update(struct.pack("<Q", len(entries)))
    for e in entries:
        h.update(struct.pack("<II", e["fg"], e["at"]))
        idxs = e["idxs"]
        h.update(struct.pack("<I", len(idxs)))
        for idx in idxs:
            h.update(struct.pack("<I", idx))
        _pack_point3d(h, e["center"], prec=6)
    return h.hexdigest()

EXPECTED = {
    "atoms": {
        "count": 671,
        "sha256": "d20618dea32021ec300dc1985a38388fc4db677c10fc64504f2151137ba031fb",
    },
    "rings": {
        "count": 10,
        "sha256": "1197f0a5dc0138c4fc3519f52b7188745de5581bec35da8ff3ebddf46c8f4f41",
    },
    "groups": {
        "count": 35,
        "sha256": "2e2823cfafdf00f1d857c09b58ac0c8255c7a2a026a34d2202bafa60aee7b8ae",
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
