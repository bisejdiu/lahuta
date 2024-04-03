"""Module for file processors."""
import functools
import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Callable, Generic, Iterable, Optional, TypeVar

from lahuta import Luni, NeighborPairs

__all__ = ["FileProcessor", "ConsumableFileProcessor"]

T = TypeVar("T")


class Worker(Generic[T]):
    """Worker strategy for file processing.

    It takes a function that processes a file and returns a dictionary of file names and results.
    It is used by the file processor classes to process files in parallel.

    Attributes:
        func (Callable[[str], T]): A function that processes a file and returns a dictionary of file names and results.

    Methods:
        execute: Execute the worker strategy.
    """

    def __init__(self, func: Callable[..., T]) -> None:
        self.func = func

    def execute(self, file_paths: list[str]) -> dict[str, T]:
        """Execute the worker strategy.

        Args:
            file_paths (list[str]): list of file paths.

        Returns:
            dict[str, str]: Dictionary of sequences.
        """
        sequence_dict = {}
        for file_path in file_paths:
            file_name = os.path.basename(file_path)
            sequence_dict[file_name] = self.func(file_path)
        return sequence_dict


class FileProcessor:
    """Class to process files and cache the results.

    Attributes:
        directory (Optional[str]): Path to a directory.
        file_list (Optional[list[str]]): list of file paths.
        allowed_file_extensions (Optional[list[str]]): list of allowed file extensions.
        n_jobs (int): Number of workers.

    Methods:
        process: Walk through the directory and process the files.

    """

    def __init__(
        self,
        directory: Optional[str] = None,
        file_list: Optional[list[str]] = None,
        allowed_file_extensions: Optional[list[str]] = None,
    ) -> None:
        self.directory = directory
        self.file_list = file_list
        self.allowed_file_extensions = set(allowed_file_extensions) if allowed_file_extensions else None

    def _is_valid_extension(self, file_name: str) -> bool:
        return not self.allowed_file_extensions or any(file_name.endswith(ext) for ext in self.allowed_file_extensions)

    def process(self, worker_func: Callable[..., T], n_jobs: int = 1, **kwargs: str) -> dict[str, T]:
        """Walk through the directory and process the files.

        Args:
            n_jobs (int): Number of workers.
            worker_func (Callable[..., T]): A function that processes a file and returns a dictionary of file names
                and results.
            *args (str): Arguments to pass to the worker function.
            **kwargs (str): Keyword arguments to pass to the worker function.

        Returns:
            dict[str, T]: Dictionary of results.

        """
        worker = Worker(functools.partial(worker_func, **kwargs))
        all_files = self._collect_file_paths()
        num_files = len(all_files)

        if num_files == 0:
            logging.warning("No files found. Nothing to do.")
            return {}

        chunk_size = max(1, num_files // n_jobs)

        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            futures = [
                executor.submit(worker.execute, all_files[i : i + chunk_size]) for i in range(0, num_files, chunk_size)
            ]

        result: dict[str, T] = {}
        for future in as_completed(futures):
            result.update(future.result())

        return result

    def _collect_file_paths(self) -> list[str]:
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

    def __repr__(self) -> str:
        files_to_show = 3
        short_file_list = (
            self.file_list[:files_to_show] + ["..."]
            if self.file_list and len(self.file_list) > files_to_show
            else self.file_list
        )
        return f"<FileProcessor(directory={self.directory}, file_list={short_file_list})>"

    def __str__(self) -> str:
        return self.__repr__()


class ConsumableFileProcessor:
    """Class to process files and cache the results.

    Attributes:
        directory (Optional[str]): Path to a directory.
        file_list (Optional[list[str]]): list of file paths.
        num_workers (int): Number of workers.
        operation (str): Operation to perform on the NeighborPairs objects.
        allowed_file_extensions (Optional[list[str]]): list of allowed file extensions.
        file_count (int): Number of files to process.

    Methods:
        walk: Walk through the directory and process the files.
        custom: Walk through the directory and process the files using a custom function.
    """

    def __init__(
        self,
        directory: Optional[str] = None,
        file_list: Optional[list[str]] = None,
        allowed_file_extensions: Optional[Iterable[str]] = None,
    ):
        self.directory = directory
        self.file_list = file_list
        self.allowed_file_extensions = allowed_file_extensions
        self.file_count = self._count_files(limit=10000)

    def _is_valid_extension(self, file_name: str) -> bool:
        return not self.allowed_file_extensions or any(file_name.endswith(ext) for ext in self.allowed_file_extensions)

    def _count_files(self, limit: int) -> int:
        count = 0
        if self.directory:
            for _, _, files in os.walk(self.directory):
                for name in files:
                    if not self._is_valid_extension(name):
                        continue
                    count += 1
                    if count >= limit:
                        break
        elif self.file_list:
            count = min(len(self.file_list), limit)
        return count

    def _worker_function(
        self,
        file_paths: list[str],
        operation: Callable[[NeighborPairs, NeighborPairs], NeighborPairs],
    ) -> NeighborPairs:
        file_iter = iter(file_paths)
        first_file_path = next(file_iter)
        luni = Luni(first_file_path)
        ns = luni.compute_neighbors()

        for file_path in file_iter:
            luni = Luni(file_path)
            ns_ = luni.compute_neighbors()
            ns = operation(ns, ns_)
        return ns

    def process(
        self, operation: Callable[[NeighborPairs, NeighborPairs], NeighborPairs], n_jobs: int = 1
    ) -> NeighborPairs | None:
        """Walk through the directory and process the files.

        Args:
            operation (Callable[[NeighborPairs, NeighborPairs], NeighborPairs]): Operation to perform on the
                NeighborPairs objects.
            n_jobs (int): Number of workers.

        Returns:
            NeighborPairs | None: NeighborPairs object.
        """
        all_files = self._collect_file_paths()
        num_files = len(all_files)

        if num_files == 0:
            return None

        chunk_size = max(1, num_files // n_jobs)
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            futures = [
                executor.submit(self._worker_function, all_files[i : i + chunk_size], operation)
                for i in range(0, num_files, chunk_size)
            ]

        futures_iter = iter(futures)
        ns = next(futures_iter).result()

        for future in futures_iter:
            ns_ = future.result()
            ns = operation(ns, ns_)

        return ns

    def _collect_file_paths(self) -> list[str]:
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

    def __repr__(self) -> str:
        return f"<FileProcessor(num_files={self.file_count})>"

    def __str__(self) -> str:
        return self.__repr__()
