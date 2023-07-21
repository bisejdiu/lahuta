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

    def __iter__(self) -> Any:
        ...

    def __new__(cls, *args: Any, **kwargs: Any) -> "AtomGroupType":
        ...

    def __getitem__(self, index: Any) -> "AtomGroupType":
        ...
