import logging
import os

from typing import Optional
from pathlib import Path
from urllib.request import urlretrieve

logging.basicConfig(level=logging.INFO)


class BaseFile:
    URL = "https://files.rcsb.org/download/"
    FILE_NAME = ""

    def __init__(self, pdb: bool = False, dir_loc: Optional[Path] = None):
        dir_loc = dir_loc or Path.cwd()
        if pdb:
            self.FILE_NAME += '.pdb'
        else:
            self.FILE_NAME += '.cif'

        self.local_path = dir_loc / self.FILE_NAME.lower()
        self.file_path = self._get_or_download()

    def _get_or_download(self) -> Path:
        if not os.path.exists(self.local_path):
            self._download_file()
        return self.local_path

    @property
    def file_loc(self) -> Path:
        return self.file_path

    def _download_file(self) -> None:
        logging.info("Downloading %s from %s", self.FILE_NAME, self.URL)
        urlretrieve(self.URL + self.FILE_NAME, self.local_path)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(file_loc={self.file_loc})"

    def __str__(self) -> str:
        return str(self.file_loc)
