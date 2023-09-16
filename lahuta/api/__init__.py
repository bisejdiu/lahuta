"""API module for the lahuta package.

It provides a high-level interface, which provides a set of classes and functions to interact with the lahuta package.
They should make it easy to use the package in a variety of contexts 
(e.g. processing files, computing neighbor pairs, etc.)
"""

from lahuta.api.np import difference, intersection, symmetric_difference, union
from lahuta.api.processors import CachedFileProcessor, FileProcessor
from lahuta.api.utils import count_unique_pairs_across_keys, download_structures, map_unique_pairs_to_keys
from lahuta.api.workers import Worker

__all__ = [
    "download_structures",
    "map_unique_pairs_to_keys",
    "count_unique_pairs_across_keys",
    "CachedFileProcessor",
    "FileProcessor",
    "Worker",
    "union",
    "intersection",
    "difference",
    "symmetric_difference",
]
