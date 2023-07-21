from typing import Literal, Optional, Protocol, Union, overload

from lahuta.core._loaders import BaseLoader
from lahuta.types.mdanalysis import AtomGroupType
from lahuta.types.openbabel import MolType


class LuniType(Protocol):
    @overload
    def to(self, fmt: Literal["mda"]) -> AtomGroupType:
        ...

    @overload
    def to(self, fmt: Literal["mol"]) -> MolType:
        ...

    def to(self, fmt: str) -> Union[MolType, AtomGroupType]:
        ...

    def to_mda(self) -> AtomGroupType:
        ...

    def to_mol(self) -> MolType:
        ...
