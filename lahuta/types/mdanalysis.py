from typing import Any, Dict, Protocol

import numpy as np
from numpy.typing import NDArray


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

    @property
    def indices(self) -> NDArray[np.int_]:
        ...

    @property
    def positions(self) -> NDArray[np.float_]:
        ...

    @property
    def universe(self) -> "UniverseType":
        ...

    @property
    def n_atoms(self) -> int:
        ...

    @property
    def names(self) -> NDArray[np.str_]:
        ...

    @property
    def ids(self) -> NDArray[np.int_]:
        ...

    @property
    def elements(self) -> NDArray[np.str_]:
        ...

    @property
    def types(self) -> NDArray[np.str_]:
        ...

    @property
    def resnames(self) -> NDArray[np.str_]:
        ...

    @property
    def resids(self) -> NDArray[np.int_]:
        ...

    @property
    def chainIDs(self) -> NDArray[np.str_]:
        ...

    def __iter__(self) -> Any:
        ...

    def __new__(cls, *args: Any, **kwargs: Any) -> "AtomGroupType":
        ...

    def __getitem__(self, index: Any) -> "AtomGroupType":
        ...

    def __len__(self) -> int:
        ...


class UniverseType(Protocol):
    @property
    def atoms(self) -> "AtomGroupType":
        ...

    @property
    def universe(self) -> "UniverseType":
        ...

    def __iter__(self) -> Any:
        ...

    def __new__(cls, *args: Any, **kwargs: Any) -> "UniverseType":
        ...
