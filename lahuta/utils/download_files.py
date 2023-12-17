"""Base class for facilitating working with PDB files."""
import logging
from pathlib import Path
from typing import Optional
from urllib.request import urlretrieve

logging.basicConfig(level=logging.INFO)


class BaseFile:
    """Base class for facilitating working with PDB files."""

    URL = "https://files.rcsb.org/download/"

    def __init__(self, pdb: bool = False, dir_loc: Optional[Path] = None):
        self.dir_loc = dir_loc or Path.cwd()
        self.file_extension = ".pdb" if pdb else ".cif"
        self.file_name = self._get_file_name()
        self.local_path = self._generate_local_path()
        self.file_path = self._get_or_download()

    def _get_file_name(self) -> str:
        """Return the file name for the instance. To be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement this method")

    def _generate_local_path(self) -> Path:
        """Generate the local path for the file."""
        return self.dir_loc / (self.file_name.lower() + self.file_extension)

    def _get_or_download(self) -> Path:
        """Get the file locally or download it if not present."""
        if not self.local_path.exists():
            self._download_file()
        return self.local_path

    def _download_file(self) -> None:
        """Download the file."""
        logging.info("Downloading %s from %s", self.file_name, self.URL)
        urlretrieve(self.URL + self.file_name + self.file_extension, self.local_path)

    @property
    def file_loc(self) -> str:
        """Return the file location."""
        return str(self.file_path)


class X2(BaseFile):
    """X2."""

    def __init__(self, pdb: bool = True) -> None:
        super().__init__(pdb=pdb)

    def _get_file_name(self) -> str:
        return "1KX2"


class Rhodopsin(BaseFile):
    """Rhodopsin."""

    def _get_file_name(self) -> str:
        return "1GZM"


class DNABound(BaseFile):
    """DNABound."""

    def _get_file_name(self) -> str:
        return "3Q2Y"
