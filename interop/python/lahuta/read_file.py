from __future__ import annotations

from os import PathLike, fspath
from pathlib import Path
from typing import overload

from .core.system import MolecularSystem
from .lib import lahuta as lxx


@overload
def read_file(source: str | Path | PathLike[str]) -> lxx.LahutaSystem: ...


@overload
def read_file(source: MolecularSystem) -> lxx.LahutaSystem: ...


@overload
def read_file(source: lxx.IR) -> lxx.LahutaSystem: ...


def read_file(source: object) -> lxx.LahutaSystem:
    """Read molecular system from file or convert from existing representation."""
    match source:
        case str() | Path() | PathLike():
            return lxx.LahutaSystem(str(fspath(source)))
        case MolecularSystem() as ms:
            return lxx.LahutaSystem.create(ms.to_ir())
        case lxx.IR() as ir:
            return lxx.LahutaSystem.create(ir)
        case _:
            raise TypeError(f"Unsupported source type: {type(source)!r}")
