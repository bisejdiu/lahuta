# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     async def f():
#         return "besian" + "sejdiu" + "@gmail.com"
#     print(asyncio.run(f()))
#
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping, Sequence

from lahuta.lib import lahuta as _lib

Database = _lib.db.Database


class DirectorySource(_lib.pipeline.sources.DirectorySource):
    def __init__(
        self,
        path: str | Path,
        *,
        recursive: bool = True,
        extensions: Iterable[str | Path] | str | Path | None = None,
        batch: int = 200,
    ) -> None:
        if extensions is None:
            ext_arg: object = ""
        elif isinstance(extensions, (str, Path)):
            ext_arg = [str(extensions)]
        else:
            ext_arg = list(dict.fromkeys(str(ext) for ext in extensions))
        super().__init__(str(path), ext_arg, bool(recursive), int(batch))


class FileSource(_lib.pipeline.sources.FileSource):
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


class LmdbSource(_lib.pipeline.sources.DatabaseSource):
    def __init__(
        self,
        db: str | Path | Database,
        *,
        database: str = "",
        batch: int = 1024,
    ) -> None:
        if isinstance(db, (str, Path)):
            super().__init__(str(db), str(database), int(batch))
        else:
            # Delegate to DatabaseHandleSource binding when a handle is provided
            _lib.pipeline.sources.DatabaseHandleSource.__init__(self, db, str(database), int(batch))


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


PipelineSource = (
    DirectorySource
    | FileSource
    | FileListSource
    | DatabaseSource
    | DatabaseHandleSource
    | LmdbSource
    | NmrSource
    | MdTrajectoriesSource
)


__all__ = [
    "PipelineSource",
    "DirectorySource",
    "FileSource",
    "FileListSource",
    "DatabaseSource",
    "DatabaseHandleSource",
    "LmdbSource",
    "NmrSource",
    "MdTrajectoriesSource",
]
