"""API module for the lahuta package.

It provides a high-level interface, which provides a set of classes and functions to interact with the lahuta package.
They should make it easy to use the package in a variety of contexts 
(e.g. processing files, computing neighbor pairs, etc.)
"""

from lahuta.api.processors import CachedFileProcessor, FileProcessor
from lahuta.api.utils import download_structures
from lahuta.api.workers import EmptyWorker, NeighborPairsComputer, SequenceComputer, WorkerStrategy

__all__ = [
    "download_structures",
    "CachedFileProcessor",
    "FileProcessor",
    "WorkerStrategy",
    "EmptyWorker",
    "NeighborPairsComputer",
    "SequenceComputer",
]
