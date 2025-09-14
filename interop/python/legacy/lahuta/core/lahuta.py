from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass
from enum import Enum, auto
from typing import Any, Protocol, TypeAlias, TypeVar, cast, overload

import numpy as np
from numpy.typing import NDArray
from typing_extensions import Self

# import numpy as np
# from numpy.typing import NDArray
from lahuta.lib._lahuta import (
    Atom,
    AtomEntityCollection,
    AtomType,
    DistanceComputation_,
    FastNS_,
    LahutaSystem,
    # MatrixType,
    NSResults_,
    Point3D,
    Residue,
    RingEntityCollection,
)

# fmt: off

IntArray:     TypeAlias = NDArray[np.int32]
StrArray:     TypeAlias = NDArray[np.str_]
FloatArray:   TypeAlias = NDArray[np.float64]
AtomPosition: TypeAlias = NDArray[np.float64] # shape=(3,)
Coordinates:  TypeAlias = NDArray[np.float64] # shape=(n_atoms, 3)

_T = TypeVar("_T")
class Residues(Protocol):
    @property
    def residues(self) -> list[Residue]: """System residues object."""; ...

    def filter(self, func: Callable[[Residue], bool]) -> Self:     """Filter residues w/ a callable."""; ...
    def map(self,    func: Callable[[Residue], _T])   -> list[_T]: """Map residues w/ a callable.""";    ...

    def __iter__(self) -> Iterator[Residue]: ...
    def __getitem__(self, index: int) -> Residue: ...


class Topology(Protocol):
    @property
    def residues(self)   -> Residues:       """System residues."""; ...
    @property
    def rings(self)      -> RingEntityCollection:    """Ring perception result."""; ...
    @property
    def atom_types(self) -> AtomEntityCollection: """Atom type assignment result."""; ...


class NeighborSearch:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, coords: list[list[float]], scale_factor: float = 1.1) -> None: ...
    @overload
    def __init__(self, coords: Coordinates,       scale_factor: float = 1.1) -> None: ...

    def __init__(self, coords: list[list[float]] | Coordinates | None = None , scale_factor: float | None = 1.1) -> None:
        """Fast Neighbor Search object."""
        self.grid = FastNS_(coords, scale_factor)
        self.is_built = False

    def build(self, cutoff: float) -> bool:
        """Build the grid."""
        self.is_built = True
        return self.grid.build(cutoff)

    def update(self, cutoff: float) -> bool:
        """Update the grid."""
        return self.grid.update(cutoff)

    def self_search(self) -> NSResults_:
        """Search for neighbors."""
        if not self.is_built:
            raise RuntimeError("Grid not built.")
        return self.grid.self_search()

    def search(self, search_coords: list[Point3D]) -> NSResults_:
        """Search for neighbors."""
        if not self.is_built:
            raise RuntimeError("Grid not built.")
        return self.grid.search(search_coords)

    def get_cutoff(self) -> float:
        """Get the cutoff."""
        return self.grid.get_cutoff()

    @staticmethod
    def dist_sq(a: Sequence[float], b: Sequence[float]) -> float:
        """Calculate the squared distance between two points."""
        return FastNS_.dist_sq(a, b)


# TODO: needs to be defined in a separate (type) file
_F_co = TypeVar("_F_co", bound=np.floating[Any], covariant=True)
MatrixType: TypeAlias = NDArray[_F_co]  # shape=(n_atoms, 3)

INT32_MAX      = np.iinfo(np.int32).max
MEM_WARN_BYTES = 4e9
MEM_WARN_SIZE  = 5_000 * 5_000

class DistanceComputation:
    @overload
    @staticmethod
    def distance(a: MatrixType[_F_co]) -> MatrixType[_F_co]: ...

    @overload
    @staticmethod
    def distance(a: Sequence[Sequence[float]]) -> MatrixType[_F_co]: ...

    @overload
    @staticmethod
    def distance(a: Sequence[float], b: Sequence[float]) -> float: ...

    @overload
    @staticmethod
    def distance(a: MatrixType[_F_co], b: MatrixType[_F_co]) -> MatrixType[_F_co]: ...

    @overload
    @staticmethod
    def distance(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]]) -> MatrixType[_F_co]: ...

    @staticmethod
    def distance(a: Sequence[float] | MatrixType[_F_co] | Sequence[Sequence[float]],
                 b: Sequence[float] | MatrixType[_F_co] | Sequence[Sequence[float]] | None = None) -> float | MatrixType[_F_co]:
        """Calculate the distance between two points or two sets of points."""
        size = len(a) * len(b) if b is not None else len(a)

        if size > INT32_MAX:
            raise ValueError(f"Cannot allocate an array with {size} elements, exceeding the max allowed {INT32_MAX} elements.")

        if size > MEM_WARN_SIZE:
            import warnings

            a_type = a.dtype if isinstance(a, np.ndarray) else type(a[0])
            b_type = a_type if b is None else b.dtype if isinstance(b, np.ndarray) else type(b[0])

            result_dtype = np.result_type(a_type, b_type)
            element_size = np.dtype(result_dtype).itemsize
            total_bytes = size * element_size

            if total_bytes > MEM_WARN_BYTES:
                warnings.warn(
                    f"Memory warning: pre-allocating an array with {size} elements "
                    f"({total_bytes / 1e9:.2f} GB) using data type {result_dtype}.",
                    stacklevel=2
                )

        # input_args = (cast(MatrixType[_F_co], a), cast(MatrixType[_F_co], b)) if b is not None else (cast(MatrixType[_F_co], a),)
        # return DistanceComputation_.distance(*input_args)
        # input_args = (a, b) if b is not None else (a,)
        # return DistanceComputation_.distance(*input_args) # type: ignore (reason): there's no point in trying to narrow down the type here

        input_args = (a, b) if b is not None else (a,)
        return DistanceComputation_.distance(*(cast(MatrixType[_F_co], x) for x in input_args)) # hack to avoid runtime overhead


    @overload
    @staticmethod
    def search(a: MatrixType[_F_co], *, cutoff: float) -> NSResults_: ...

    @overload
    @staticmethod
    def search(a: Sequence[Sequence[float]], *, cutoff: float) -> NSResults_: ...

    @overload
    @staticmethod
    def search(a: MatrixType[_F_co], b: MatrixType[_F_co], *, cutoff: float) -> NSResults_: ...

    @overload
    @staticmethod
    def search(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]], *, cutoff: float) -> NSResults_: ...

    @staticmethod
    def search(a: MatrixType[_F_co] | Sequence[Sequence[float]],
               b: MatrixType[_F_co] | Sequence[Sequence[float]] | float | None = None,
               *,
               cutoff: float = 5.0) -> NSResults_:
        """Search for neighbors."""
        if b is None:
            return DistanceComputation_.search(a, cutoff)
        # a, b = (np.array(x) if not isinstance(x, np.ndarray) else x for x in (a, b)) # works, but creats a new array

        # bit of a hack, but we avoid any runtime overhead.
        return DistanceComputation_.search(cast(MatrixType[_F_co], a), cast(MatrixType[_F_co], b), cutoff)


class Lahuta:
    def __init__(self, file_name: str) -> None:
        self.file_name = file_name
        self.luni = LahutaSystem(file_name)

    @classmethod
    def from_file(cls, file_name: str) -> "Lahuta": """Create a Lahuta object from a file."""; return cls(file_name)

    def get_conformer_positions(self, conformer_index: int = -1) -> Coordinates:
        """Get atom positions for a conformer."""
        return self.luni.get_positions(conformer_index)

    def get_topology(self) -> Topology:
        """Get the topology of the system."""
        return self.luni.get_topology()


    def find_neighbors(self, radius: float, res_diff: int) -> NSResults_:
        return self.luni.find_neighbors(radius, res_diff)


    @property
    def indices(self)    -> IntArray:     """Atom indices.""";    return self.luni.indices
    @property
    def atom_nums(self)  -> IntArray:     """Atomic numbers.""";  return self.luni.atom_nums
    @property
    def resids(self)     -> IntArray:     """Residue IDs.""";     return self.luni.resids
    @property
    def resindices(self) -> IntArray:     """Residue indices."""; return self.luni.resindices
    @property
    def names(self)      -> StrArray:     """Atom names.""";      return self.luni.names
    @property
    def symbols(self)    -> StrArray:     """Atom symbols.""";    return self.luni.symbols
    @property
    def elements(self)   -> StrArray:     """Element symbols."""; return self.luni.elements
    @property
    def resnames(self)   -> StrArray:     """Residue names.""";   return self.luni.resnames
    @property
    def chainlabels(self)-> IntArray:     """Chain labels.""";    return self.luni.chainlabels
    @property
    def positions(self)  -> AtomPosition: """Atom positions.""";  return self.luni.positions # FIX: this doesn't return AtomPosition but rather Coordinates
    @property
    def n_atoms(self)    -> int:          """Number of atoms."""; return self.luni.n_atoms

    def get_atom(self, index: int) -> Atom: """Get an atom by index."""; return self.luni.get_atom(index)
