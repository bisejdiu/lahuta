from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping, Sequence

from lahuta.lib import lahuta as _lib

Database = _lib.db.Database


class Source(_lib.pipeline.sources.Source):
    """Base class for all pipeline sources."""


class DirectorySource(_lib.pipeline.sources.DirectorySource):
    def __init__(self, path: str | Path, *, ext: str = "", recursive: bool = True, batch: int = 200) -> None:
        super().__init__(str(path), str(ext), bool(recursive), int(batch))


class FilesSource(_lib.pipeline.sources.FilesSource):
    def __init__(self, files: Sequence[str | Path] | str | Path) -> None:
        if isinstance(files, (str, Path)):
            payload = [str(files)]
        else:
            payload = [str(f) for f in files]
        super().__init__(payload)


class FileListSource(_lib.pipeline.sources.FileListSource):
    def __init__(self, path: str | Path) -> None:
        super().__init__(str(path))


class DatabaseSource(_lib.pipeline.sources.DatabaseSource):
    def __init__(self, path: str | Path, *, database: str = "", batch: int = 1024) -> None:
        super().__init__(str(path), str(database), int(batch))


class DatabaseHandleSource(_lib.pipeline.sources.DatabaseHandleSource):
    def __init__(self, db: Database, *, database: str = "", batch: int = 1024) -> None:
        super().__init__(db, str(database), int(batch))


class NmrSource(_lib.pipeline.sources.NmrSource):
    def __init__(self, files: Sequence[str | Path] | str | Path) -> None:
        if isinstance(files, (str, Path)):
            payload = [str(files)]
        else:
            payload = [str(f) for f in files]
        super().__init__(payload)


class MdTrajectoriesSource(_lib.pipeline.sources.MdTrajectoriesSource):
    def __init__(self, trajectories: Iterable[Mapping[str, object] | Sequence[object]]) -> None:
        super().__init__(trajectories)


__all__ = [
    "Source",
    "DirectorySource",
    "FilesSource",
    "FileListSource",
    "DatabaseSource",
    "DatabaseHandleSource",
    "NmrSource",
    "MdTrajectoriesSource",
]
