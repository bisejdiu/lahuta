"""Worker strategies for the API."""
import os
from typing import Callable, Generic, TypeVar

__all__ = ["Worker"]

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

    def __init__(self, func: Callable[[str], T]) -> None:
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
