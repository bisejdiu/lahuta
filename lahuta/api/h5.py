"""Module for handling HDF5 files."""
try:
    import h5py
except ImportError:
    raise ImportError("h5py is not installed. Please install it to use this module.") from None

import tempfile
from multiprocessing import Lock
from pathlib import Path

import numpy as np
from numpy.typing import NDArray


class HDF5Handler:
    """Class for handling HDF5 files.

    Args:
        file_location (str, optional): Location of the HDF5 file. If not provided, a temporary file will be created.
        dtype (str | NDArray[np.void], optional): Data type of the data to be saved. Defaults to None.

    Attributes:
        file_location (str): Location of the HDF5 file.
        dtype (str | NDArray[np.void]): Data type of the data to be saved.
    """

    def __init__(self, file_location: str | Path | None = None, dtype: str | NDArray[np.void] | None = None):
        if file_location:
            self.file_location = Path(file_location)
        else:
            temp_file = tempfile.NamedTemporaryFile(delete=True)
            self.file_location = Path(temp_file.name)

        self.dtype = dtype or [
            ("chainids", "<U25"),
            ("names", "<U25"),
            ("resids", "<U25"),
            ("resnames", "<U25"),
            ("gns", "<U25"),
        ]
        self.lock = Lock()

    def _convert_dtype(self, data: NDArray[np.void]) -> NDArray[np.void]:
        # Convert the string fields to byte strings for HDF5 compatibility
        if isinstance(self.dtype, str):
            return np.array(data, dtype=self.dtype)

        new_dtype = []
        for name, dtype in self.dtype:
            new_dtype.append((name, "S25" if dtype == "<U25" else dtype))

        return np.array(data, dtype=new_dtype)

    def save_data(self, data: dict[str, NDArray[np.void]]) -> None:
        """Save data to the HDF5 file.

        Args:
            data (dict[str, NDArray[np.void]]): Data to be saved.
        """
        with h5py.File(self.file_location, "w") as h5file:
            for key, value in data.items():
                h5file.create_dataset(key, data=self._convert_dtype(value))

    def append_data(self, key: str, value: NDArray[np.void]) -> None:
        """Append data to the HDF5 file. This operation is thread-safe.

        Args:
            key (str): Key for the data to be appended.
            value (NDArray[np.void]): Data to be appended.
        """
        with self.lock, h5py.File(self.file_location, "a") as h5file:
            if key in h5file:
                raise KeyError(
                    f"Dataset with key '{key}' already exists. Use a different key or delete the existing dataset."
                )
            h5file.create_dataset(key, data=self._convert_dtype(value))

    def read_data(self, key: str) -> NDArray[np.void]:
        """Read data from the HDF5 file with the given key.

        Args:
            key (str): Key for the data to be read.

        Raises:
            KeyError: If the key is not found in the HDF5 file.

        Returns:
            NDArray[np.void]: Data read from the HDF5 file.
        """
        with h5py.File(self.file_location, "r") as h5file:
            if key not in h5file:
                raise KeyError(f"Key {key} not found in HDF5 file.")
            encoded_array = h5file[key][:]
            return np.array(encoded_array, self.dtype)

    def load_data(self) -> dict[str, NDArray[np.void]]:
        """Load data from the HDF5 file.

        Depending on the size of the file, this may take a while and use a lot of memory.

        Returns:
            dict[str, NDArray[np.void]]: Data loaded from the HDF5 file.
        """
        data = {}
        with h5py.File(self.file_location, "r") as h5file:
            for key in h5file:
                encoded_array = h5file[key][:]
                data[key] = np.array(encoded_array, self.dtype)
        return data
