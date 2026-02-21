# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     print(functools.reduce(lambda a, b: a + b, ["besian", "sejdiu", "@gmail.com"], ""))
#
from __future__ import annotations

from pathlib import Path
from typing import TypeAlias

from .lib import lahuta as _lib
from .pipeline.wrapper import Pipeline
from .sources import DatabaseHandleSource, DirectorySource

Database: TypeAlias = _lib.db.Database


class LahutaDB:
    """
    - Provides creation from a directory of AlphaFold models
    - Exposes read operations
    - Offers a convenience method to create a Pipeline supported by this DB
    """

    def __init__(self, path: str | Path):
        self._path = str(path)
        self._db = _lib.db.Database(self._path)

    @classmethod
    def from_handle(cls, path: str | Path, db: Database) -> LahutaDB:
        """Wrap an existing LMDB handle without reopening the environment."""
        obj = cls.__new__(cls)
        obj._path = str(path)
        obj._db = db
        return obj

    # TODO: too arg-heavy
    @classmethod
    def create_from_directory(
        cls,
        directory: str | Path,
        out: str | Path,
        *,
        ext: str = ".cif.gz",
        recursive: bool = False,
        batch: int = 256,
        threads: int = 4,
    ) -> "LahutaDB":
        """Create a model database from a directory of inputs."""
        source = DirectorySource(directory, extensions=[str(ext)] if str(ext) else None, recursive=bool(recursive), batch=int(batch))
        mgr = _lib.pipeline.StageManager(source)

        # Emit serialized ModelRecord on channel 'db'
        pack = _lib.pipeline.ModelPackTask("db")
        mgr.add_task("pack", [], pack, True)

        db = _lib.db.Database(str(out))
        sink = _lib.pipeline.LmdbSink(db, int(batch))
        mgr.connect_sink("db", sink)

        mgr.set_auto_builtins(True)
        mgr.run(int(threads))

        # Reuse the same handle to avoid reopening the environment immediately
        return cls.from_handle(out, db)

    def keys(self) -> list[str]:
        return self._db.keys()

    def get_model(self, key: str) -> _lib.LahutaSystem:
        return self._db.get_model(key)

    def to_pipeline(self, *, batch: int = 256) -> Pipeline:
        """Create a Pipeline reading from this DB with model mode auto configured."""
        source = DatabaseHandleSource(self._db, batch=int(batch))
        p = Pipeline(source)
        return p

    @property
    def path(self) -> str:
        return self._path


__all__ = ["LahutaDB"]
