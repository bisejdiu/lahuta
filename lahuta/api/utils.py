"""Utility functions for the API."""
from collections import defaultdict
from enum import Enum
from operator import itemgetter
from pathlib import Path
from typing import Any, Literal, Optional

import numpy as np

from lahuta.tests.base import BaseFile

__all__ = ["download_structures", "count_unique_pairs_across_keys", "map_unique_pairs_to_keys"]


class URLs(Enum):
    """Enum class for URLs."""

    RCSB = "https://files.rcsb.org/download/"
    AlphaFold = "https://alphafold.ebi.ac.uk/files/"


def download_structures(
    pdb_ids: list[str],
    url: str | URLs = URLs.RCSB,
    pdb_or_cif: Literal["pdb", "cif"] = "pdb",
    dir_loc: Optional[str | Path] = None,
) -> dict[str, str]:
    """Download PDB/CIF files from the RCSB PDB database.

    Args:
        pdb_ids (list[str]): List of PDB IDs to download.
        url (URLs, optional): URL to download from. Defaults to URLs.RCSB.
        pdb_or_cif (Literal["pdb", "cif"], optional): File format to download. Defaults to "pdb".
        dir_loc (Path): Path to directory where PDB files will be downloaded.

    Returns:
        dict[str, str]: Dictionary with PDB IDs as keys and file locations as values.
    """
    dir_loc = Path.cwd() if dir_loc is None else Path(dir_loc)
    url_value = url.value if isinstance(url, URLs) else url

    if not dir_loc.exists():
        dir_loc.mkdir(parents=True)

    pdb_file_locations: dict[str, str] = {}
    for pdb_id in pdb_ids:
        tf_class = type("PDBDownloader_" + pdb_id, (BaseFile,), {"FILE_NAME": pdb_id, "URL": url_value})
        tf = tf_class(pdb=pdb_or_cif == "pdb", dir_loc=dir_loc)
        pdb_file_locations[pdb_id] = tf.file_loc.as_posix()

    return pdb_file_locations


def count_unique_pairs_across_keys(
    data_dict: dict[str, Any], mapping: Optional[dict[str, str]] = None
) -> dict[str, int]:
    """Count the number of unique pairs across all keys in a dictionary.

    Args:
        data_dict (dict[str, Any]): Dictionary with keys as keys and values as iterables of pairs.
        mapping (dict[str, str], optional): Mapping to apply to the pairs. Defaults to None.

    Returns:
        dict[str, int]: Dictionary with unique pairs as keys and their counts as values.
    """
    result: defaultdict[str, int] = defaultdict(int)

    for pairs in data_dict.values():
        if mapping:
            pairs = np.vectorize(mapping.get)(pairs)  # noqa: PLW2901
        unique_pairs = {tuple(i) for i in pairs}

        for pair in unique_pairs:
            result[str(pair)] += 1

    # Sort by value, highest first
    return dict(sorted(result.items(), key=itemgetter(1), reverse=True))


def map_unique_pairs_to_keys(
    data_dict: dict[str, Any], mapping: Optional[dict[str, str]] = None
) -> dict[str, list[str]]:
    """Map unique pairs to the keys in which they occur.

    Args:
        data_dict (dict[str, Any]): Dictionary with keys as keys and values as iterables of pairs.
        mapping (dict[str, str], optional): Mapping to apply to the pairs. Defaults to None.

    Returns:
        dict[str, list[str]]: Dictionary with unique pairs as keys and lists of keys as values.
    """
    result = defaultdict(list)

    for key, pairs in data_dict.items():
        if mapping:
            pairs = np.vectorize(mapping.get)(pairs)  # noqa: PLW2901
        unique_pairs = {tuple(i) for i in pairs}

        for pair in unique_pairs:
            result[str(pair)].append(key)

    # Sort by the length of the value list, highest first
    return dict(sorted(result.items(), key=lambda x: len(x[1]), reverse=True))
