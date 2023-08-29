"""Module for file processors."""
import logging
import os
import sys
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Generic, List, Optional, Type, TypeVar, cast

from lahuta import Luni, NeighborPairs

T = TypeVar("T")


class WorkerStrategy(ABC, Generic[T]):
    """Abstract base class for worker strategies.

    Methods:
        execute: Execute the worker strategy.
    """

    @abstractmethod
    def execute(self, file_paths: List[str]) -> dict[str, T]:
        """Execute the worker strategy.

        Args:
            file_paths (List[str]): List of file paths.

        Returns:
            dict[str, T]: Dictionary of objects.
        """
        ...


class NeighborPairsComputer(WorkerStrategy[NeighborPairs]):
    """Worker strategy to compute neighbor pairs.

    Methods:
        execute: Execute the worker strategy.
    """

    def execute(self, file_paths: List[str]) -> dict[str, NeighborPairs]:
        """Execute the worker strategy.

        Args:
            file_paths (List[str]): List of file paths.

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

    def execute(self, file_paths: List[str]) -> dict[str, str]:
        """Execute the worker strategy.

        Args:
            file_paths (List[str]): List of file paths.

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

    def execute(self, file_paths: List[str]) -> dict[str, None]:
        """Execute the worker strategy.

        Args:
            file_paths (List[str]): List of file paths.

        Returns:
            dict[str, None]: Dictionary of None values.
        """
        return {file_name: None for file_name in file_paths}


class CachedFileProcessor(Generic[T]):
    """Class to process files and cache the results.

    Attributes:
        directory (Optional[str]): Path to a directory.
        file_list (Optional[List[str]]): List of file paths.
        allowed_file_extensions (Optional[List[str]]): List of allowed file extensions.
        num_workers (int): Number of workers.
        worker (WorkerStrategy[T]): A worker strategy.

    Methods:
        walk: Walk through the directory and process the files.
        results: Get the cached results.
    """

    def __init__(
        self,
        directory: Optional[str] = None,
        file_list: Optional[List[str]] = None,
        allowed_file_extensions: Optional[List[str]] = None,
        num_workers: int = 4,
        worker: Type[WorkerStrategy[T]] = cast(Type[WorkerStrategy[T]], NeighborPairsComputer),  # noqa: B008
    ):
        self.worker = worker()
        self.num_workers = num_workers
        self.directory = directory
        self.file_list = file_list
        self.allowed_file_extensions = set(allowed_file_extensions) if allowed_file_extensions else None
        self._cache: dict[str, Any] = {}

    def _is_valid_extension(self, file_name: str) -> bool:
        return not self.allowed_file_extensions or any(file_name.endswith(ext) for ext in self.allowed_file_extensions)

    def process(self) -> None:
        """Walk through the directory and process the files."""
        all_files = self._collect_file_paths()
        num_files = len(all_files)

        if num_files == 0:
            logging.warning("No files found. Nothing to do.")
            return

        chunk_size = max(1, num_files // self.num_workers)

        with ThreadPoolExecutor(max_workers=self.num_workers) as executor:
            futures = [
                executor.submit(self.worker.execute, all_files[i : i + chunk_size])
                for i in range(0, num_files, chunk_size)
            ]

        for future in futures:
            self._cache.update(future.result())

    def _collect_file_paths(self) -> List[str]:
        if self.directory:
            return [
                os.path.join(root, name)
                for root, _, files in os.walk(self.directory)
                for name in files
                if self._is_valid_extension(name)
            ]
        if self.file_list:
            return [file for file in self.file_list if self._is_valid_extension(file)]

        raise ValueError("Either directory or file_list must be provided")

    @property
    def results(self) -> dict[str, T]:
        """Get the cached results."""
        if not self._cache:
            logging.warning("No results available. Did you forget to call walk()?")
        return self._cache

    def __repr__(self) -> str:
        num_objects = len(self._cache)
        memory_footprint = sum(sys.getsizeof(ns) for ns in self._cache.values())
        return f"<CachedFileProcessor(num_objects={num_objects}, memory_footprint={memory_footprint} bytes)>"

    def __str__(self) -> str:
        return self.__repr__()
