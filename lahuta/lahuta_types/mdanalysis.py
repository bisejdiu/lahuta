"""
Module: mdanalysis.py

This module provides indirect typing support for MDAnalysis, a library for 
the analysis of molecular dynamics simulations. Despite the absence of inherent 
typing support in MDAnalysis, the `mdanalysis.py` module uses Protocols to define 
interfaces that allow type hinting for MDAnalysis objects.

While currently essential, this module is projected to become obsolete as soon 
as the MDAnalysis library incorporates typing natively, eliminating the need for 
indirect typing.

Classes:
    ResidueGroupType(Protocol): Provides a typing interface for MDAnalysis ResidueGroup objects.
    AtomGroupType(Protocol): Provides a typing interface for MDAnalysis AtomGroup objects.
    UniverseType(Protocol): Provides a typing interface for MDAnalysis Universe objects.

Example:
    def compute_distance(group: AtomGroupType) -> float:
        # An example function that operates on an mda.AtomGroup object.
        # The AtomGroupType Protocol allows type hinting for the group parameter.
        ...

    universe: UniverseType = ...
    atoms = universe.atoms
    compute_distance(atoms)
"""

from typing import Any, Optional, Protocol, Union

import numpy as np
from numpy.typing import NDArray


# pylint: disable=missing-function-docstring
class ResidueGroupType(Protocol):
    """A typing interface for MDAnalysis ResidueGroup objects."""

    @property
    def atoms(self) -> "AtomGroupType":
        ...

    @property
    def resnames(self) -> NDArray[np.str_]:
        ...

    def __iter__(self) -> Any:
        ...


# pylint: disable=missing-function-docstring
class AtomGroupType(Protocol):
    """A typing interface for MDAnalysis AtomGroup objects."""

    def select_atoms(self, selection: str) -> "AtomGroupType":
        ...

    @property
    def residues(self) -> "ResidueGroupType":
        ...

    @property
    def atoms(self) -> "AtomGroupType":
        ...

    @property
    def indices(self) -> NDArray[np.int32]:
        ...

    @property
    def positions(self) -> NDArray[np.float32]:
        ...

    @positions.setter
    def positions(self, positions: NDArray[np.float32]) -> None:
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
    def ids(self) -> NDArray[np.int32]:
        ...

    @property
    def ix(self) -> NDArray[np.int32]:
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
    def resname(self) -> NDArray[np.str_]:
        ...

    @property
    def resids(self) -> NDArray[np.int32]:
        ...

    @property
    def resindices(self) -> NDArray[np.int_]:
        ...

    @property
    def chainIDs(self) -> NDArray[np.str_]:  # pylint: disable=invalid-name
        ...

    @property
    def vdw_radii(self) -> NDArray[np.float16]:
        ...

    def copy(self) -> "AtomGroupType":
        ...

    def __iter__(self) -> Any:
        ...

    def __new__(cls, *args: Any, **kwargs: Any) -> "AtomGroupType":
        ...

    def __getitem__(self, index: Union[int, slice, NDArray[np.int32]]) -> "AtomGroupType":
        ...

    def __len__(self) -> int:
        ...


# pylint: disable=missing-function-docstring
class UniverseType(Protocol):
    """A typing interface for MDAnalysis Universe objects."""

    @property
    def atoms(self) -> "AtomGroupType":
        ...

    @property
    def universe(self) -> "UniverseType":
        ...

    @property
    def dimensions(self) -> Optional[NDArray[np.float32]]:
        ...

    @property
    def trajectory(self) -> "TrajectoryType":
        ...

    def select_atoms(self, selection: str) -> "AtomGroupType":
        ...

    def copy(self) -> "UniverseType":
        ...

    def add_TopologyAttr(self, name: str, attr: Any) -> None:  # pylint: disable=invalid-name
        ...

    def __iter__(self) -> Any:
        ...

    def __new__(cls, *args: Any, **kwargs: Any) -> "UniverseType":
        ...


class TrajectoryType(Protocol):
    @property
    def n_frames(self) -> int:
        ...

    def __getitem__(self, index: Union[int, slice, NDArray[np.int32]]) -> "TimeStepType":
        ...

    def __iter__(self) -> "TimeStepType":
        ...


class TimeStepType(Protocol):
    @property
    def frame(self) -> int:
        ...

    def __iter__(self) -> "TimeStepType":
        ...

    def __next__(self) -> "TimeStepType":
        ...
