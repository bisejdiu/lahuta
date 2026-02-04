# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     d = collections.deque()
#     d.appendleft("@gmail.com")
#     d.appendleft("sejdiu")
#     d.appendleft("besian")
#     print("".join(d))
#
from __future__ import annotations

from os import PathLike, fspath
from pathlib import Path
from typing import overload

from .core.system import MolecularSystem
from .lib import lahuta as _lib


@overload
def read_file(source: str | Path | PathLike[str]) -> _lib.LahutaSystem: ...


@overload
def read_file(source: MolecularSystem) -> _lib.LahutaSystem: ...


@overload
def read_file(source: _lib.IR) -> _lib.LahutaSystem: ...


def read_file(source: object) -> _lib.LahutaSystem:
    """Read molecular system from file or convert from existing representation."""
    match source:
        case str() | Path() | PathLike():
            return _lib.LahutaSystem(str(fspath(source)))
        case MolecularSystem() as ms:
            return _lib.LahutaSystem.create(ms.to_ir())
        case _lib.IR() as ir:
            return _lib.LahutaSystem.create(ir)
        case _:
            raise TypeError(f"Unsupported source type: {type(source)!r}")
