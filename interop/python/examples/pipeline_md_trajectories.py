from __future__ import annotations

import os
from pathlib import Path

from lahuta import ContactProvider, InteractionType
from lahuta.pipeline import (
    ContactTask,
    FileOutput,
    InMemoryPolicy,
    OutputFormat,
    Pipeline,
    PipelineContext,
)
from lahuta.sources import MdTrajectoriesSource

DATA_DIR = Path(__file__).resolve().parents[3] / "old_files_and_dirs" / "gpcrmd"
STRUCTURE_FILE = DATA_DIR / "10828_dyn_85.pdb"
TRAJECTORY_FILE = DATA_DIR / "10824_trj_85.xtc"


def md_contacts_basic() -> None:
    """
    Basic MD trajectory processing with contacts computation.
    Uses tuple format: (structure_path, [trajectory_paths])
    """

    # tuple format: (structure_path, [xtc_path1, xtc_path2, ...])
    trajectories = [
        (str(STRUCTURE_FILE), [str(TRAJECTORY_FILE)]),
    ]

    source = MdTrajectoriesSource(trajectories)
    p = Pipeline(source)

    # Contacts computation will be run for each frame in the trajectory
    p.add_task(
        name="contacts",
        task=ContactTask(provider=ContactProvider.MolStar, interaction_type=InteractionType.All),
        in_memory_policy=InMemoryPolicy.Keep,
    )

    out = p.run(threads=os.cpu_count() or 8)

    #
    # Can also use out.to_dict("contacts"), but it's much slower as many Python objects are created
    # with few optimization opportunities.  - Besian, November 26, 2025
    #
    contacts = out.to_numpy("contacts")
    if contacts:
        first = contacts[0]
        frame_index = first.get("frame_index", "N/A")
        n_contacts = first.get("num_contacts", "N/A")
        print(f"First frame index: {frame_index}")
        print(f"Number of contacts in first frame: {n_contacts}")
        print(f"Available data keys: {list(first.keys())}")
    print()


def md_contacts_dict_format() -> None:
    """
    MD trajectory processing using dict format for trajectory specification.
    Dict keys: 'structure' (or 'gro'), 'xtc' (or 'xtcs'), optional 'id'/'session_id'
    """

    # dict format with explicit keys
    trajectories = [
        {
            "structure": str(STRUCTURE_FILE),
            "xtc": str(TRAJECTORY_FILE),
            "session_id": "gpcrmd_85",  # Optional custom identifier
        },
    ]

    source = MdTrajectoriesSource(trajectories)
    p = Pipeline(source)

    p.add_task(
        name="contacts",
        task=ContactTask(
            provider=ContactProvider.MolStar,
            interaction_type=InteractionType.HydrogenBond,  # Only H-bonds
        ),
        in_memory_policy=InMemoryPolicy.Keep,
    )

    out = p.run(threads=os.cpu_count() or 8)
    contacts = out.to_dict("contacts")
    print(f"Processed {len(contacts)} frames (H-bonds only)")
    print()


def md_contacts_to_file() -> None:
    """MD trajectory processing with output to JSONL file."""

    out_path = Path("md_contacts.jsonl")
    if out_path.exists():
        out_path.unlink()

    trajectories = [(str(STRUCTURE_FILE), [str(TRAJECTORY_FILE)])]
    source = MdTrajectoriesSource(trajectories)
    p = Pipeline(source)

    p.add_task(
        name="contacts",
        task=ContactTask(),  # molstar provider, all interactions by default
        out=[FileOutput(out_path, fmt=OutputFormat.JSON)],
        in_memory_policy=InMemoryPolicy.Drop,  # Don't keep in memory
    )

    p.run(threads=os.cpu_count() or 8)
    files = p.file_outputs()
    print(f"Output files: {files.get('contacts')}")
    if out_path.exists():
        print(f"File size: {out_path.stat().st_size / (1024 * 1024):.2f} MB")
    print()


def md_contacts_with_analysis() -> None:
    """MD trajectory processing with custom Python analysis task."""

    trajectories = [(str(STRUCTURE_FILE), [str(TRAJECTORY_FILE)])]
    source = MdTrajectoriesSource(trajectories)
    p = Pipeline(source)

    p.add_task(
        name="contacts",
        task=ContactTask(),
        in_memory_policy=InMemoryPolicy.Keep,
    )

    # Custom analysis that depends on contacts
    def analyze_frame(ctx: PipelineContext) -> dict:
        """Summarize contact information for each frame."""
        return {
            "frame_id": ctx.path,
            "processed": True,
        }

    p.add_task(
        name="frame_summary",
        task=analyze_frame,
        depends=["contacts"],
        in_memory_policy=InMemoryPolicy.Keep,
        thread_safe=True,  # it's a pure function, so safe for concurrent execution
    )

    out = p.run(threads=os.cpu_count() or 8)

    contacts = out.to_numpy("contacts")
    summaries = out.to_dict("frame_summary")
    print(f"Contacts: {len(contacts)} frames")
    print(f"Summaries: {len(summaries)} frames")
    print()


def md_interaction_types() -> None:
    """Filter specific interaction types from MD trajectory."""

    trajectories = [(str(STRUCTURE_FILE), [str(TRAJECTORY_FILE)])]
    source = MdTrajectoriesSource(trajectories)
    p = Pipeline(source)

    # Hydrogen bonds only
    p.add_task(
        name="hbonds",
        task=ContactTask(interaction_type=InteractionType.HydrogenBond),
        channel="hbonds",
        in_memory_policy=InMemoryPolicy.Keep,
    )

    # Ionic interactions only
    p.add_task(
        name="ionic",
        task=ContactTask(interaction_type=InteractionType.Ionic),
        channel="ionic",
        in_memory_policy=InMemoryPolicy.Keep,
    )

    out = p.run(threads=os.cpu_count() or 8)
    for channel in ["hbonds", "ionic"]:
        results = out.to_dict(channel)
        print(f"{channel}: {len(results)} frames")
    print()


if __name__ == "__main__":
    if not STRUCTURE_FILE.exists() or not TRAJECTORY_FILE.exists():
        print("Required data files are missing.")
        exit(1)

    md_contacts_basic()
    md_contacts_dict_format()
    md_contacts_to_file()
    md_contacts_with_analysis()
    md_interaction_types()
