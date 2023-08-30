"""Utility functions for the API."""
from pathlib import Path
from typing import Literal, Optional

from lahuta.tests.base import BaseFile

__all__ = ["download_structures"]


def download_structures(
    pdb_ids: list[str], pdb_or_cif: Literal["pdb", "cif"] = "pdb", dir_loc: Optional[str | Path] = None
) -> dict[str, str]:
    """Download PDB/CIF files from the RCSB PDB database.

    Args:
        pdb_ids (list[str]): List of PDB IDs to download.
        pdb_or_cif (Literal["pdb", "cif"], optional): File format to download. Defaults to "pdb".
        dir_loc (Path): Path to directory where PDB files will be downloaded.

    Returns:
        dict[str, str]: Dictionary with PDB IDs as keys and file locations as values.
    """
    dir_loc = Path.cwd() if dir_loc is None else Path(dir_loc)

    if not dir_loc.exists():
        dir_loc.mkdir(parents=True)

    pdb_file_locations: dict[str, str] = {}
    for pdb_id in pdb_ids:
        tf_class = type("PDBDownloader_" + pdb_id, (BaseFile,), {"FILE_NAME": pdb_id})
        tf = tf_class(pdb=pdb_or_cif == "pdb", dir_loc=dir_loc)
        pdb_file_locations[pdb_id] = tf.file_loc.as_posix()

    return pdb_file_locations
