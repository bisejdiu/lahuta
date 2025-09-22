import os
from collections import Counter
from pathlib import Path
from typing import Any, Callable

import pytest

from lahuta import AtomTypingMethod, TopologyBuildingOptions
from lahuta.db import LahutaDB
from lahuta.lib import lahuta as lxx
from lahuta.pipeline import InMemoryPolicy
from lahuta.pipeline.tasks import ContactTask
from lahuta.pipeline.wrapper import Pipeline

# Expected contact counts for 1kx2_small.cif using Arpeggio
EXPECTED_ARPEGGIO_CONTACT_COUNTS: dict[str, int] = {
    "SulphurPi": 6,
    "CarbonPi": 10,
    "DonorPi": 2,
    "CationPi": 0,
    "Carbonyl": 8,
    "Ionic": 11,
    "Aromatic": 20,
    "Hydrophobic": 150,
    "VanDerWaals": 56,
    "HydrogenBond": 64,
    "WeakHydrogenBond": 44,
    "PolarHydrogenBond": 100,
    "WeakPolarHydrogenBond": 93,
}


def _validate_contact_counts(
    contacts: list[Any], expected_counts: dict[str, int], contact_type_extractor: Callable[[Any], str]
) -> None:
    """Validate contact counts against expected values."""

    type_counts = Counter(contact_type_extractor(c) for c in contacts)

    for contact_type, expected_count in expected_counts.items():
        actual_count = type_counts.get(contact_type, 0)
        assert actual_count == expected_count, (
            f"Contact type '{contact_type}': expected {expected_count}, got {actual_count}"
        )

    total_expected = sum(expected_counts.values())
    total_actual = len(contacts)
    assert total_actual == total_expected, f"Total contacts: expected {total_expected}, got {total_actual}"


@pytest.mark.parametrize("provider", [lxx.ContactProvider.MolStar, lxx.ContactProvider.Arpeggio])
def test_create_db_and_compute_contacts(tmp_path: Path, provider: lxx.ContactProvider) -> None:
    """Test creating a database from directory and computing contacts using different providers."""
    data_dir = Path("core/data/models")
    assert data_dir.exists(), "Test models directory missing"

    db_path = tmp_path / "models_lmdb"
    dbw = LahutaDB.create_from_directory(data_dir, db_path, ext=".cif.gz", recursive=False, batch=50, threads=2)
    db = dbw._db

    keys = db.keys()
    assert len(keys) == 2

    p_read = Pipeline.from_database_handle(db)

    task_name = "contacts_molstar" if provider == lxx.ContactProvider.MolStar else "contacts_arpeggio"
    p_read.add_task(name=task_name, task=ContactTask(provider=provider))

    out = p_read.run(threads=2)
    assert task_name in out
    recs = out[task_name]
    assert len(recs) == 2

    if provider == lxx.ContactProvider.MolStar:
        expected = {
            "AF-P0CL56-F1-model_v4.cif.gz": 57,
            "AF-Q57552-F1-model_v4.cif.gz": 687,
        }
    else:
        expected = {
            "AF-P0CL56-F1-model_v4.cif.gz": 186,
            "AF-Q57552-F1-model_v4.cif.gz": 2243,
        }

    got: dict[str, int] = {}
    for rec in recs:
        assert isinstance(rec, dict)
        assert "file_path" in rec and "num_contacts" in rec and "provider" in rec
        base = os.path.basename(rec["file_path"])
        got[base] = int(rec["num_contacts"])

    assert got == expected


def test_db_pipeline_exposes_model_mode(tmp_path: Path) -> None:
    data_dir = Path("core/data/models")
    assert data_dir.exists(), "Test models directory missing"

    db_path = tmp_path / "models_lmdb"
    dbw = LahutaDB.create_from_directory(data_dir, db_path, ext=".cif.gz", recursive=False, batch=50, threads=2)
    db = dbw._db

    p_read = Pipeline.from_database_handle(db)

    def get_mode(ctx: Any) -> dict[str, str]:
        sys = ctx.get_system()
        return {"mode": "model" if getattr(sys, "is_model", False) else "structure"}

    p_read.add_task(name="mode", task=get_mode, in_memory_policy=InMemoryPolicy.Keep)
    out = p_read.run(threads=2)
    modes = [rec["mode"] for rec in out["mode"]]
    assert modes == ["model", "model"]


def test_arpeggio_contacts_file_pipeline() -> None:
    data_file = Path("core/data/1kx2_small.cif")
    assert data_file.exists(), "Test file missing"

    p = Pipeline.from_files(str(data_file))
    p.add_task(name="contacts", task=ContactTask(provider=lxx.ContactProvider.Arpeggio))

    out = p.run(threads=1)
    assert "contacts" in out
    recs = out["contacts"]
    assert len(recs) == 1

    rec = recs[0]
    assert isinstance(rec, dict)
    assert "file_path" in rec and "contacts" in rec and "provider" in rec
    assert str(rec["provider"]).lower() == "arpeggio"

    contacts = rec.get("contacts", [])
    _validate_contact_counts(contacts, EXPECTED_ARPEGGIO_CONTACT_COUNTS, lambda c: str(c.get("type", "")))


def test_arpeggio_contacts_lahuta_system() -> None:
    from lahuta.entities import compute_contacts

    data_file = Path("core/data/1kx2_small.cif")
    assert data_file.exists(), "Test file missing"

    opts = TopologyBuildingOptions()
    opts.atom_typing_method = AtomTypingMethod.Arpeggio
    sys = lxx.LahutaSystem(str(data_file))
    if not sys.build_topology(opts):
        raise RuntimeError("Failed to build topology")

    top = sys.get_topology()

    contacts = compute_contacts(top, provider="arpeggio")

    contacts_list = list(contacts)
    _validate_contact_counts(contacts_list, EXPECTED_ARPEGGIO_CONTACT_COUNTS, lambda c: str(c.type))
