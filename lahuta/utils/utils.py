"""Utility functions for the lahuta package."""
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import IO, Generator


def string_to_file(text: str, extension: str) -> str:
    """Create a temporary file with the specified extension.

    Args:
        text (str): The text to write to the file.
        extension (str): The extension of the file.

    Returns:
        str: The path to the temporary file.
    """
    # Create a temporary file with the specified extension
    temp_file = tempfile.NamedTemporaryFile(mode="w+t", suffix="." + extension, delete=False)
    temp_file_path = temp_file.name
    temp_file.write(text)
    temp_file.close()

    return temp_file_path


@contextmanager
def managed_temp_file(text: str, extension: str, delete: bool = True) -> Generator[IO[str], None, None]:
    """Context manager that creates a temporary file with the given text and extension.
    The file is deleted after exiting the context if the delete parameter is True.
    """
    temp_file = tempfile.NamedTemporaryFile(mode="w+t", suffix="." + extension, delete=False)
    temp_file_path = Path(temp_file.name)
    try:
        temp_file.write(text)
        temp_file.flush()
        temp_file.seek(0)
        yield temp_file
    finally:
        # Close the file after the context is exited
        temp_file.close()
        # Delete the file if the delete flag is True
        if delete and temp_file_path.exists():
            temp_file_path.unlink()
