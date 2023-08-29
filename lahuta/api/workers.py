"""Worker strategies for the API."""
import os
from abc import ABC, abstractmethod
from typing import Generic, TypeVar

from lahuta.core.luni import Luni
from lahuta.core.neighbors import NeighborPairs

__all__ = ["WorkerStrategy", "NeighborPairsComputer", "SequenceComputer", "EmptyWorker"]

T = TypeVar("T")


class WorkerStrategy(ABC, Generic[T]):
    """Abstract base class for worker strategies.

    Methods:
        execute: Execute the worker strategy.
    """

    @abstractmethod
    def execute(self, file_paths: list[str]) -> dict[str, T]:
        """Execute the worker strategy.

        Args:
            file_paths (list[str]): list of file paths.

        Returns:
            dict[str, T]: Dictionary of objects.
        """
        ...


class NeighborPairsComputer(WorkerStrategy[NeighborPairs]):
    """Worker strategy to compute neighbor pairs.

    Methods:
        execute: Execute the worker strategy.
    """

    def execute(self, file_paths: list[str]) -> dict[str, NeighborPairs]:
        """Execute the worker strategy.

        Args:
            file_paths (list[str]): list of file paths.

        Returns:
            dict[str, NeighborPairs]: Dictionary of NeighborPairs objects.
        """
        ns_dict = {}
        for file_path in file_paths:
            file_name = os.path.basename(file_path)
            luni = Luni(file_path)
            ns = luni.compute_neighbors()
            ns_dict[file_name] = ns
        return ns_dict


class SequenceComputer(WorkerStrategy[str]):
    """Worker strategy to compute sequences.

    Methods:
        execute: Execute the worker strategy.
    """

    def execute(self, file_paths: list[str]) -> dict[str, str]:
        """Execute the worker strategy.

        Args:
            file_paths (list[str]): list of file paths.

        Returns:
            dict[str, str]: Dictionary of sequences.
        """
        sequence_dict = {}
        for file_path in file_paths:
            file_name = os.path.basename(file_path)
            luni = Luni(file_path)
            sequence_dict[file_name] = luni.sequence
        return sequence_dict


class EmptyWorker(WorkerStrategy[None]):
    """Worker strategy to do nothing.

    Methods:
        execute: Execute the worker strategy.
    """

    def execute(self, file_paths: list[str]) -> dict[str, None]:
        """Execute the worker strategy.

        Args:
            file_paths (list[str]): list of file paths.

        Returns:
            dict[str, None]: Dictionary of None values.
        """
        return {file_name: None for file_name in file_paths}
