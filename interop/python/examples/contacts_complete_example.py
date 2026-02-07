# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     print(f"{'besian'}{'sejdiu'}@{'gmail.com'}")
#
"""
The output format matches what `lahuta contacts -f <file>` produces:
{
    "file_path": "...",
    "success": true,
    "provider": "molstar",
    "contact_type": "All",
    "num_contacts": N,
    "frame_index": 0,
    "contacts": [
        {"lhs": "...", "rhs": "...", "distance_sq": ..., "type": "..."},
        ...
    ]
}

Entity Format:
- Atoms:  "atomIdx-atomName-resSeq-resName-chainId" e.g., "4-O-1-MET-A"
- Groups/Rings: "(atomIdx-atomName, atomIdx-atomName, ...)-resSeq-resName-chainId"
         e.g., "(128-NE, 129-NH1, 130-NH2)-14-ARG-A"

This example relies on Contact.to_dict(topology), which uses the same compact
entity formatting as the CLI.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, TypedDict

from lahuta import (
    InteractionType,
    LahutaSystem,
)
from lahuta.entities import compute_contacts
from lahuta.entities.api import ContactProviderType
from lahuta.pipeline import InMemoryPolicy, Pipeline, PipelineContext
from lahuta.sources import DirectorySource

# fmt: off
DATA_DIR = Path(__file__).resolve().parents[3] / "core" / "data"
DATA = DATA_DIR / "models" / "AF-P0CL56-F1-model_v4.cif.gz"

class ContactDict(TypedDict):
    lhs: str
    rhs: str
    distance_sq: float
    type: str


class ContactOutputDict(TypedDict):
    file_path: str
    success: bool
    provider: str
    contact_type: str
    num_contacts: int
    frame_index: int
    contacts: list[ContactDict]


def compute_and_resolve_contacts(
    file_path: str | Path,
    provider: ContactProviderType = "molstar",
    only: InteractionType | None = None
) -> ContactOutputDict:
    """Compute and resolve contacts."""
    file_path = Path(file_path)

    system = LahutaSystem(str(file_path))
    assert system.build_topology(), "Failed to build topology"

    top = system.get_topology()
    # compute contacts and resolve entities
    contact_set = compute_contacts(top, provider=provider, only=only)
    contacts_list: list[ContactDict] = [contact.to_dict(top) for contact in contact_set]

    return {
        "file_path": str(file_path),
        "success": True,
        "provider": provider,
        "contact_type": "All" if only is None else str(only),
        "num_contacts": len(contacts_list),
        "frame_index": 0,
        "contacts": contacts_list,
    }

def print_contact_summary(result: ContactOutputDict, max_contacts: int = 10) -> None:
    """Print a summary of the computed contacts."""
    print(f"File: {result['file_path']}")
    print(f"Success: {result['success']}")
    print(f"Provider: {result['provider']}")
    print(f"Contact Type: {result['contact_type']}")
    print(f"Total Contacts: {result['num_contacts']}")
    print()

    if not result.get("contacts"):
        print("No contacts found.")
        return

    print(f"First {min(max_contacts, len(result['contacts']))} contacts:")
    for contact in result["contacts"][:max_contacts]:
        print(f"  {contact['lhs']:20} -&- {contact['rhs']:20} {contact['distance_sq']:8.3f} A^2 ({contact['type']})")


def pipeline_basic_directory(
    directory: Path | None = None,
    provider: ContactProviderType = "molstar",
    only: InteractionType | None = None,
) -> list[ContactOutputDict]:
    """
    Process a directory of structure files using the Pipeline, computing
    contacts on the Python side (bypassing C++ ContactTask).
    """

    def compute_contacts_task(ctx: PipelineContext) -> dict[str, Any]:
        top = ctx.get_topology()
        contact_set = compute_contacts(top, provider=provider, only=only)
        contacts_list: list[ContactDict] = [contact.to_dict(top) for contact in contact_set]

        return {
            "file_path": ctx.path,
            "success": True,
            "provider": provider,
            "contact_type": "All" if only is None else str(only),
            "num_contacts": len(contacts_list),
            "frame_index": 0,
            "contacts": contacts_list,
        }

    source = DirectorySource(directory or DATA_DIR, recursive=False, extensions=[".cif"], batch=64)
    p = Pipeline(source)
    p.add_task(
        name="contacts",
        task=compute_contacts_task,
        in_memory_policy=InMemoryPolicy.Keep,
    )

    out = p.run(threads=4)
    return out.to_dict("contacts")


def pipeline_custom_postprocessing(directory: Path | None = None, provider: ContactProviderType = "molstar") -> list[dict[str, Any]]:
    """
    Process a directory with custom post-processing: group contacts by type
    and compute summary statistics (count, min/max/avg distance_sq per type).
    """

    def contacts_with_stats(ctx: PipelineContext) -> dict[str, Any]:
        top = ctx.get_topology()
        contact_set = compute_contacts(top, provider=provider, only=None)

        by_type: dict[str, list[dict]] = {}
        for contact in contact_set:
            record = contact.to_dict(top)
            ctype = record["type"]
            if ctype not in by_type:
                by_type[ctype] = []
            by_type[ctype].append(
                {
                    "lhs": record["lhs"],
                    "rhs": record["rhs"],
                    "distance_sq": record["distance_sq"],
                }
            )

        stats = {}
        for ctype, contacts in by_type.items():
            distances_sq = [c["distance_sq"] for c in contacts]
            stats[ctype] = {
                "count": len(contacts),
                # these checks are very defensive, as distances should never be empty
                "min_distance_sq": min(distances_sq) if distances_sq else 0,
                "max_distance_sq": max(distances_sq) if distances_sq else 0,
                "avg_distance_sq": sum(distances_sq) / len(distances_sq) if distances_sq else 0,
            }

        return {
            "file_path": ctx.path,
            "success": True,
            "total_contacts": sum(len(v) for v in by_type.values()),
            "stats_by_type": stats,
            "contacts_by_type": by_type,
        }

    source = DirectorySource(directory or DATA_DIR, recursive=False, extensions=[".cif"], batch=64)
    p = Pipeline(source)
    p.add_task(
        name="contacts_stats",
        task=contacts_with_stats,
        in_memory_policy=InMemoryPolicy.Keep,
    )

    out = p.run(threads=4)
    return out.to_dict("contacts_stats")


def main() -> None:
    if not DATA.exists():
        print(f"Data file not found: {DATA}")
        print("Please update the DATA path to point to a valid structure file.")
        return

    print("\nComputing all contacts...")
    result = compute_and_resolve_contacts(DATA)
    print_contact_summary(result)

    print("JSON output (first contact only):")
    result_preview = result.copy()
    result_preview["contacts"] = result["contacts"][:1] if result["contacts"] else []
    print(json.dumps(result_preview, indent=2))
    print()

    if result["contacts"]:
        from collections import Counter

        type_counts = Counter(c["type"] for c in result["contacts"])
        print("Contact type distribution:")
        for ctype, count in sorted(type_counts.items(), key=lambda x: -x[1]):
            print(f"  {ctype}: {count}")


    print("Pipeline with basic directory processing")
    results = pipeline_basic_directory()
    print(f"\nProcessed {len(results)} files:")
    for r in results[:5]:
        print(f"  {Path(r['file_path']).name}: {r['num_contacts']} contacts")
    if len(results) > 5:
        print(f"  ... and {len(results) - 5} more")


    print("Pipeline with custom post-processing with stats")
    stats_results = pipeline_custom_postprocessing()
    print(f"\nProcessed {len(stats_results)} files with stats:")
    for r in stats_results[:3]:
        print(f"\n  {Path(r['file_path']).name}: {r['total_contacts']} total contacts")
        print("    Stats by type:")
        for ctype, stats in r.get("stats_by_type", {}).items():
            print(f"      {ctype}: {stats['count']} (avg: {stats['avg_distance_sq']:.2f} A^2)")
    if len(stats_results) > 3:
        print(f"\n  ... and {len(stats_results) - 3} more")


if __name__ == "__main__":
    main()
