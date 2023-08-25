import logging
import os
from pathlib import Path
from urllib.request import urlretrieve


def find_project_root(start_from: str = __file__, marker: str = 'pyproject.toml', max_depth: int = 3) -> Path:
    path = Path(start_from).resolve().parent
    depth = 0
    while depth < max_depth:
        if (path / marker).exists():
            return path
        path = path.parent
        depth += 1
    raise FileNotFoundError(f"Project root with marker '{marker}' not found within {max_depth} levels.")


class BaseFile:
    URL = "https://files.rcsb.org/download/"
    FILE_NAME = ""

    def __init__(self, pdb: bool = False):
        project_root = find_project_root() / 'tests' / 'data'
        if pdb:
            self.FILE_NAME += '.pdb'
        else:
            self.FILE_NAME += '.cif'

        self.local_path = project_root / self.FILE_NAME.lower()
        self.file_path = self._get_or_download()

    def _get_or_download(self) -> Path:
        if not os.path.exists(self.local_path):
            print('downloading')
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
