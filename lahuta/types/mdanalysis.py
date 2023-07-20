from typing import Any, Dict, Protocol


class ResidueGroupType(Protocol):
    @property
    def atoms(self) -> "AtomGroupType":
        ...

    def __iter__(self) -> Any:
        ...


class AtomGroupType(Protocol):
    def select_atoms(self, selection: str) -> "AtomGroupType":
        ...

    @property
    def residues(self) -> "ResidueGroupType":
        ...

    @property
    def atoms(self) -> "AtomGroupType":
        ...

    def __iter__(self) -> Any:
        ...
