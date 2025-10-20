from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

import pytest

import lahuta as lxx
from lahuta.pipeline import Pipeline
from lahuta.pipeline.tasks import ContactTask
from lahuta.sources import FileSource


# fmt: off
def _load_expected_contacts(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _compute_getcontacts_contacts(structure_path: Path) -> dict[str, Any]:
    pipeline = Pipeline(FileSource(str(structure_path)))
    pipeline.add_task(
        name="gc",
        task=ContactTask(provider=lxx.ContactProvider.GetContacts),
    )
    output = pipeline.run(threads=1)
    records = output.to_dict("gc")

    assert isinstance(records, list), "Expected records to be a list"
    assert len(records) == 1, f"Expected exactly one record, got {len(records)}"

    return records[0]


def _count_by_interaction_type(contacts_data: dict[str, Any]) -> Counter[str]:
    contacts = contacts_data.get("contacts", [])
    return Counter(str(contact.get("type", "")) for contact in contacts)


def _validate_basic_metadata(computed: dict[str, Any], expected: dict[str, Any], structure_name: str) -> None:
    assert computed.get("success") is True, f"{structure_name}: Computation failed"

    assert str(computed.get("provider", "")).lower() == "getcontacts", (
        f"{structure_name}: Expected provider 'getcontacts', got '{computed.get('provider')}'"
    )

    assert computed.get("provider")    == expected.get("provider"), f"{structure_name}: Provider mismatch"
    assert computed.get("frame_index") == expected.get("frame_index"), f"{structure_name}: Frame index mismatch"


def _validate_interaction_counts(computed_counts: Counter[str], expected_counts: Counter[str], structure_name: str) -> None:
    computed_types = set(computed_counts.keys())
    expected_types = set(expected_counts.keys())

    missing_types = expected_types - computed_types
    extra_types   = computed_types - expected_types

    assert not missing_types, f"{structure_name}: Missing interaction types: {sorted(missing_types)}"
    assert not extra_types,   f"{structure_name}: Unexpected interaction types: {sorted(extra_types)}"

    mismatches: list[str] = []
    for interaction_type in sorted(expected_types):
        computed_count = computed_counts[interaction_type]
        expected_count = expected_counts[interaction_type]

        if computed_count != expected_count:
            mismatches.append(f"  {interaction_type}: computed={computed_count}, expected={expected_count}")

    assert not mismatches, f"{structure_name}: Interaction count mismatches:\n" + "\n".join(mismatches)


def _validate_total_contact_count(computed: dict[str, Any], expected: dict[str, Any], structure_name: str) -> None:
    computed_total = len(computed.get("contacts", []))
    expected_total = len(expected.get("contacts", []))

    assert computed_total == expected_total, (
        f"{structure_name}: Total contact count mismatch: computed={computed_total}, expected={expected_total}"
    )


def test_getcontacts_contact_count_sanity_check(data_dir: Path) -> None:
    structure_path = data_dir / "1kx2_small.cif"

    if not structure_path.exists():
        pytest.skip(f"Test structure not found: {structure_path}")

    computed_data = _compute_getcontacts_contacts(structure_path)

    total_contacts = len(computed_data.get("contacts", []))
    assert total_contacts > 100, (
        f"Expected more than 100 contacts for 1kx2_small.cif with GetContacts, got {total_contacts}"
    )

    counts = _count_by_interaction_type(computed_data)
    assert len(counts) >= 3, f"Expected at least 3 interaction types from GetContacts, got {len(counts)}"

    expected_types = {"HydrogenBond", "VanDerWaals"}
    found_types = expected_types & set(counts.keys())
    assert found_types == expected_types, f"Expected interaction types {expected_types}, found {found_types}"


@pytest.mark.parametrize(
    "structure_file,expected_file",
    [
        ("1kx2_small.cif", "1kx2_getcontacts.json"),
        ("1tqn.cif", "1tqn_getcontacts.json"),
        ("models/AF-Q57552-F1-model_v4.cif.gz", "Q57552_getcontacts.json"),
    ],
)
def test_getcontacts_contacts_regression(data_dir: Path, structure_file: str, expected_file: str) -> None:
    """
    Test that GetContacts contact computation produces consistent results.

    We do not check:
    - Exact ordering of contacts
    - Exact distance values
    - Exact residue string formatting
    """
    structure_path = data_dir / structure_file
    expected_path  = data_dir / "contacts" / expected_file

    assert structure_path.exists(), f"Structure file not found: {structure_path}"
    assert expected_path.exists(),  f"Expected contacts file not found: {expected_path}"

    structure_name = structure_path.name

    expected_data = _load_expected_contacts(expected_path)
    computed_data = _compute_getcontacts_contacts(structure_path)

    _validate_basic_metadata(computed_data, expected_data, structure_name)
    _validate_total_contact_count(computed_data, expected_data, structure_name)

    computed_counts = _count_by_interaction_type(computed_data)
    expected_counts = _count_by_interaction_type(expected_data)

    _validate_interaction_counts(computed_counts, expected_counts, structure_name)
