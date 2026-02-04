# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     d = collections.deque()
#     d.appendleft("@gmail.com")
#     d.appendleft("sejdiu")
#     d.appendleft("besian")
#     print("".join(d))
#
import os
from pathlib import Path

from lahuta import ContactProvider, InteractionType
from lahuta.pipeline import (
    ContactTask,
    FileOutput,
    InMemoryPolicy,
    OutputFormat,
    Pipeline,
)
from lahuta.sources import DirectorySource

DATA_DIR = Path("/Users/bsejdiu/data/UP000000805_243232_METJA_v4/")


def compute_contacts_all_providers(output_prefix: str = "contacts", threads: int = 4) -> None:
    """
    Compute contacts using all three providers in a single pipeline run.
    Each provider's results are saved to a separate NDJSON file.
    """
    source = DirectorySource(DATA_DIR, recursive=False, extensions=[".cif.gz"], batch=64)

    p = Pipeline(source)
    p.params("system").is_model = True  # working with AF2 models

    arpeggio_file = Path(f"{output_prefix}_arpeggio.ndjson")
    molstar_file = Path(f"{output_prefix}_molstar.ndjson")
    getcontacts_file = Path(f"{output_prefix}_getcontacts.ndjson")

    # Clean up existing output files. Lahuta does not do this and will append to existing files.
    for file in [arpeggio_file, molstar_file, getcontacts_file]:
        if file.exists():
            file.unlink()

    p.add_task(
        name="arpeggio_contacts",
        task=ContactTask(provider=ContactProvider.Arpeggio),
        out=[FileOutput(arpeggio_file, fmt=OutputFormat.JSON)],
        in_memory_policy=InMemoryPolicy.Drop,
    )

    p.add_task(
        name="molstar_contacts",
        task=ContactTask(provider=ContactProvider.MolStar, interaction_type=InteractionType.All),
        out=[FileOutput(molstar_file, fmt=OutputFormat.JSON)],
        in_memory_policy=InMemoryPolicy.Drop,
    )

    p.add_task(
        name="getcontacts_contacts",
        task=ContactTask(provider=ContactProvider.GetContacts, interaction_type=InteractionType.All),
        out=[FileOutput(getcontacts_file, fmt=OutputFormat.JSON)],
        in_memory_policy=InMemoryPolicy.Drop,
    )

    print(f"Processing files from: {DATA_DIR}")
    p.run(threads=threads)

    files = p.file_outputs()
    print("\nResults saved:")
    print(f"  Arpeggio:    {files.get('arpeggio_contacts')}")
    print(f"  MolStar:     {files.get('molstar_contacts')}")
    print(f"  GetContacts: {files.get('getcontacts_contacts')}")


if __name__ == "__main__":
    compute_contacts_all_providers(output_prefix="contacts_all", threads=os.cpu_count() or 4)
