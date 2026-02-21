# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         a = "besian"
#         b = "sejdiu"
#         c = "@gmail.com"
#     print(vars(Email)["a"] + vars(Email)["b"] + vars(Email)["c"])
#
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

import pytest

from lahuta import ContactProvider
from lahuta.pipeline import Pipeline
from lahuta.pipeline.tasks import ContactTask
from lahuta.sources import FileSource


# fmt: off
def _load_expected_contacts(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _compute_molstar_contacts(structure_path: Path) -> dict[str, Any]:
    pipeline = Pipeline(FileSource(str(structure_path)))
    pipeline.add_task(
        name="ms",
        task=ContactTask(provider=ContactProvider.MolStar),
    )
    output = pipeline.run(threads=1)
    records = output.to_dict("ms")

    assert isinstance(records, list), "Expected records to be a list"
    assert len(records) == 1, f"Expected exactly one record, got {len(records)}"

    return records[0]


def _count_by_interaction_type(contacts_data: dict[str, Any]) -> Counter[str]:
    contacts = contacts_data.get("contacts", [])
    return Counter(str(contact.get("type", "")) for contact in contacts)


def _validate_basic_metadata(computed: dict[str, Any], expected: dict[str, Any], structure_name: str) -> None:
    assert computed.get("success") is True, f"{structure_name}: Computation failed"
    assert str(computed.get("provider", "")).lower() == "molstar", (
        f"{structure_name}: Expected provider 'molstar', got '{computed.get('provider')}'"
    )
    assert computed.get("provider")    == expected.get("provider"),    f"{structure_name}: Provider mismatch"
    assert computed.get("frame_index") == expected.get("frame_index"), f"{structure_name}: Frame index mismatch"


def _validate_interaction_counts(computed_counts: Counter[str], expected_counts: Counter[str], structure_name: str) -> None:
    # Check that we have the same set of interaction types
    computed_types = set(computed_counts.keys())
    expected_types = set(expected_counts.keys())

    missing_types = expected_types - computed_types
    extra_types   = computed_types - expected_types

    assert not missing_types, f"{structure_name}: Missing interaction types: {sorted(missing_types)}"
    assert not extra_types,   f"{structure_name}: Unexpected interaction types: {sorted(extra_types)}"

    # Check counts for each type
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


def test_molstar_contact_count_sanity_check(data_dir: Path) -> None:
    structure_path = data_dir / "1kx2_small.cif"

    if not structure_path.exists():
        pytest.skip(f"Test structure not found: {structure_path}")

    computed_data = _compute_molstar_contacts(structure_path)

    total_contacts = len(computed_data.get("contacts", []))
    assert total_contacts > 50, f"Expected more than 50 contacts for 1kx2_small.cif, got {total_contacts}"

    counts = _count_by_interaction_type(computed_data)

    assert len(counts) >= 3, f"Expected at least 3 interaction types, got {len(counts)}"

    assert "HydrogenBond" in counts, "Expected HydrogenBond interactions"
    assert "Hydrophobic"  in counts, "Expected Hydrophobic interactions"


@pytest.mark.parametrize(
    "structure_file, expected_file",
    [
        ("1kx2_small.cif", "1kx2_molstar.json"),
        ("1tqn.cif", "1tqn_molstar.json"),
        ("models/AF-Q57552-F1-model_v4.cif.gz", "Q57552_molstar.json"),
    ],
)
def test_molstar_contacts_regression(data_dir: Path, structure_file: str, expected_file: str) -> None:
    """
    Test that MolStar contact computation produces consistent results.

    We do not check:
    - Exact ordering of contacts
    - Exact distance values
    - Exact residue string formatting
    """
    structure_path = data_dir / structure_file
    expected_path  = data_dir / "contacts" / expected_file

    # Verify files exist
    assert structure_path.exists(), f"Structure file not found: {structure_path}"
    assert expected_path.exists(),  f"Expected contacts file not found: {expected_path}"

    structure_name = structure_path.name

    expected_data = _load_expected_contacts(expected_path)
    computed_data = _compute_molstar_contacts(structure_path)

    _validate_basic_metadata(computed_data, expected_data, structure_name)
    _validate_total_contact_count(computed_data, expected_data, structure_name)

    computed_counts = _count_by_interaction_type(computed_data)
    expected_counts = _count_by_interaction_type(expected_data)

    _validate_interaction_counts(computed_counts, expected_counts, structure_name)
