from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum, IntEnum
from typing import (
    Any,
    Callable,
    ClassVar,
    Generic,
    Iterable,
    Iterator,
    Literal,
    NewType,
    Protocol,
    Self,
    TypeAlias,
    TypeVar,
    overload,
)

import numpy as np
from numpy.typing import NDArray

from lahuta.lib._contacts import Contact as Contact
from lahuta.lib._contacts import Contacts as Contacts
from lahuta.lib._contacts import EntityID as EntityID
from lahuta.lib._contacts import EntityType as EntityType
from lahuta.lib._contacts import InteractionOptions as InteractionOptions
from lahuta.lib._contacts import Interactions as Interactions
from lahuta.lib._contacts import InteractionType as InteractionType
from lahuta.lib._entities import AtomEntity as AtomEntity
from lahuta.lib._entities import AtomEntityCollection as AtomEntityCollection
from lahuta.lib._entities import EntityNeighborSearch as EntityNeighborSearch
from lahuta.lib._entities import GroupEntity as GroupEntity
from lahuta.lib._entities import GroupEntityCollection as GroupEntityCollection
from lahuta.lib._entities import RingEntity as RingEntity
from lahuta.lib._entities import RingEntityCollection as RingEntityCollection

# from lahuta.lib import atom_types as atom_types
# from lahuta.lib.atom_types import AtomType as AtomType

# from lahuta.lib._lahuta.atom_types import AtomType

# __all__ = ["CCAtomType"]

# fmt: off

# NOTE: it would be nice to simplify the typing by encapsulating the `list`
# and `NDArray` types in a custom type alias or protocol, but type checkers
# do not expand these, so the user is left wondering what StrIter is.
# class StrIter(Protocol):
#     def __iter__(self) -> Iterator[str]: ...
# StrIter: TypeAlias = list[str] | NDArray[np.str_]

Labels:  TypeAlias = list[str] | NDArray[np.str_]
Indices: TypeAlias = list[int] | NDArray[np.int32]

Vector: TypeAlias = NDArray[np.float64] # shape=(3,)
Matrix: TypeAlias = NDArray[np.float64] # shape=(n_atoms, 3)
# Vector = NewType("Vector", NDArray[np.float64])
# Matrix = NewType("Matrix", NDArray[np.float64])

Bytes: TypeAlias = int

_Float_co = TypeVar("_Float_co", bound=np.floating[Any], covariant=True)
VectorType: TypeAlias = NDArray[_Float_co]  # shape=(3,)
MatrixType: TypeAlias = NDArray[_Float_co]  # shape=(n_atoms, 3)

class Atom:
    @property
    def atomic_number(self) -> int: ...
    @property
    def symbol(self) -> str: ...
    @property
    def idx(self) -> int: ...
    @property
    def degree(self) -> int: ...
    @property
    def formal_charge(self) -> int: ...
    @property
    def hybridization(self) -> int: ...
    @property
    def num_explicit_Hs(self) -> int: ...
    @property
    def num_implicit_Hs(self) -> int: ...
    @property
    def total_num_Hs(self) -> int: ...
    @property
    def total_valence(self) -> int: ...
    @property
    def explicit_valence(self) -> int: ...
    @property
    def implicit_valence(self) -> int: ...
    @property
    def is_aromatic(self) -> bool: ...
    @property
    def mass(self) -> float: ...


class Point3D:
    x: float
    y: float
    z: float
    def __init__(self, x: float = 0.0, y: float = 0.0, z: float = 0.0) -> None: ...



class NSResults_:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, pairs: list[tuple[int, int]], distances: list[float]) -> None: ...
    @overload
    def __init__(self, pairs: NDArray[np.int32], distances: NDArray[np.float32]) -> None: ...
    def add(self) -> None: ...
    def add_neighbors(self) -> None: ...
    def reserve_space(self) -> None: ...
    def size(self) -> int: ...
    def clear(self) -> None: ...
    def get_pairs(self) -> NDArray[np.int32]: ...
    def get_distances(self) -> NDArray[np.float32]: ...
    @overload
    def filter(self, value: float) -> NSResults_: ...
    @overload
    def filter(self, indices: list[int]) -> NSResults_: ... # replace list[int] with AtomIndices or something
    @overload
    def filter(self, indices: list[int], column_idx: int) -> NSResults_: ... # replace list[int] with AtomIndices or something


class FastNS_:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, coords: NDArray[np.float64] | None, scale_factor: float | None = 1.1) -> None: ...
    @overload
    def __init__(self, coords: list[list[float]] | None, scale_factor: float | None = 1.1) -> None: ...

    def build(self, cutoff: float) -> bool: ...
    def update(self, cutoff: float) -> bool: ...
    def adaptive_build(self, cutoff: float, max_retries: int = 50) -> bool: ...
    def self_search(self) -> NSResults_: ...
    def search(self, search_coords: list[Point3D]) -> NSResults_: ...
    def get_cutoff(self) -> float: ...
    @staticmethod
    def dist_sq(a: Sequence[float], b: Sequence[float]) -> float: ...


class DistanceComputation_:
    @overload
    @staticmethod
    def distance(a: MatrixType[_Float_co])     -> MatrixType[_Float_co]: ...

    @overload
    @staticmethod
    def distance(a: Sequence[Sequence[float]]) -> MatrixType[_Float_co]: ...

    @overload
    @staticmethod
    def distance(a: Sequence[float],           b: Sequence[float]) -> float: ...

    @overload
    @staticmethod
    def distance(a: MatrixType[_Float_co],     b: MatrixType[_Float_co])         -> MatrixType[_Float_co]: ...

    @overload
    @staticmethod
    def distance(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]]) -> MatrixType[_Float_co]: ...

    @overload
    @staticmethod
    def search(a: MatrixType[_Float_co],         cutoff: float) -> NSResults_: ...

    @overload
    @staticmethod
    def search(a: Sequence[Sequence[float]], cutoff: float) -> NSResults_: ...

    @overload
    @staticmethod
    def search(a: MatrixType[_Float_co],         b: MatrixType[_Float_co],         cutoff: float) -> NSResults_: ...

    @overload
    @staticmethod
    def search(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]], cutoff: float) -> NSResults_: ...



# Intermediate Representation
@dataclass
class IR:
    atom_indices: list[int] | NDArray[np.int32]
    atom_nums:    list[int] | NDArray[np.int32]
    resids:       list[int] | NDArray[np.int32]
    atom_names:   list[str] | NDArray[np.str_]
    resnames:     list[str] | NDArray[np.str_]
    chainlabels:  list[str] | NDArray[np.str_]
    positions:    NDArray[np.float32]

class Residue:
    chain_id: str
    number:   int
    name:     str
    alt_loc:  str
    atoms:    list[int]


_T = TypeVar("_T")
class Residues_:
    @property
    def residues(self) -> list[Residue]: ...
    def residue_map(self)  -> dict[int, list[Residue]]: ...
    def filter(self, func: Callable[[Residue], bool]) -> Residues_: ...
    def map(self,    func: Callable[[Residue], _T])   -> list[_T]: ...
    def total_size(self) -> Bytes: ...
    def __iter__(self) -> Iterator[Residue]: ...
    def __getitem__(self, index: int) -> Residue: ...


class ContactComputerType(Enum):
    None_    = 0
    Arpeggio = 1
    Molstar  = 2


@dataclass
class TopologyBuildingOptions:
    bonded_search_cutoff: float = 4.5 # the maximum distance to search for bonded atoms
    atom_typing_method: ContactComputerType = ContactComputerType.Molstar


@dataclass
class Topology_:
    @property
    def residues(self)   -> Residues_: ...
    @property
    def rings(self)      -> RingEntityCollection: ...
    @property
    def atom_types(self) -> AtomEntityCollection: ...
    def total_size(self)  -> Bytes: ...
    # TODO: indexing by integer should give an atom object
    # TODO: indexing by list | slice should give a list of atom objects


# TODO: main lahuta class should support creating an MDAnalysis object
# TODO: MDAnalysis installation should maybe be optional? (not sure)
class LahutaCPP:
    # @overload
    # def __init__(self) -> None: ...
    @overload
    def __init__(self, filename: str, contact_type: int = 1) -> None: ...
    @overload
    def __init__(self, ir: IR) -> None: ...
    def build_topology(self, t_ops: TopologyBuildingOptions | None = None) -> bool: ...
    # def get_positions(self) -> list[list[float]]: ...
    def get_positions(self, conformer_index: int = -1) -> NDArray[np.float64]: ...
    # def filter(self) -> LahutaCPP: ...
    # def filter_luni(self, indices: list[int]) -> LahutaCPP: ...
    def find_neighbors(self, radius: float, res_diff: int) -> Any: ...
    def find_ring_neighbors(self, radius: float, res_diff: int) -> Any: ...
    # def find_neighbors(self, radius: float, res_diff: int) -> AtomAtomNeighbors: ...
    # def find_ring_neighbors(self, radius: float, res_diff: int) -> AtomRingNeighbors: ...
    def get_atom_types(self) -> list[AtomType]: ...
    def get_cutoff(self) -> float: ...
    def get_rings(self) -> RingEntityCollection: ...
    # def n_atoms(self) -> int: ...
    def parse_expression(self, expression: str) -> bool: ...
    # def get_indices(self) -> NDArray[np.int32]: ...
    def get_indices(self) -> list[int]: ...
    def get_atomic_numbers(self) -> list[int]: ...
    def get_elements(self) -> list[str]: ...
    def get_symbols(self) -> list[str]: ...
    def get_names(self) -> list[str]: ...
    def get_resids(self) -> list[int]: ...
    def get_resindices(self) -> list[int]: ...
    def get_resnames(self) -> list[str]: ...
    def get_chainlabels(self) -> list[str]: ...

    def get_atom(self, index: int) -> Atom: ...

    def total_size(self) -> Bytes: ...

    def get_topology(self) -> Topology_: ...
    @staticmethod
    def factorize(labels: list[str] | NDArray[np.str_]) -> list[int] | NDArray[np.int32]: ...
    @overload
    @staticmethod
    def count_unique(labels: list[str] | NDArray[np.str_]) -> int: ...
    @overload
    @staticmethod
    def count_unique(labels: list[int] | NDArray[np.int32]) -> int: ...

    # FIX: this actually returns a list of strings
    @staticmethod
    def find_elements(atomic_numbers: list[int] | NDArray[np.int32]) -> NDArray[np.str_]: ...
    @property
    def file_name(self) -> str: ...
    @property
    def n_atoms(self) -> int: ...
    @property
    def indices(self) -> NDArray[np.int32]: ...
    @property
    def atom_nums(self) -> NDArray[np.int32]: ...
    @property
    def elements(self) -> NDArray[np.str_]: ...
    @property
    def names(self) -> NDArray[np.str_]: ...
    @property
    def resids(self) -> NDArray[np.int32]: ...
    @property
    def resindices(self) -> NDArray[np.int32]: ...
    @property
    def resnames(self) -> NDArray[np.str_]: ...
    @property
    def symbols(self) -> NDArray[np.str_]: ...
    @property
    def chainlabels(self) -> NDArray[np.int32]: ...
    @property
    def positions(self) -> NDArray[np.float64]: ...


# FIX: point is a coordinate point; the typing should reflect that
class RingData:
    # atom_ids: list[int]
    def __init__(self) -> None: ...
    # def compute_angle(self, point: list[float]) -> float: ...
    def atom_ids(self) -> list[int]: ...
    @property
    def center(self) -> NDArray[np.floating]: ...
    @property
    def norm(self) -> NDArray[np.floating]: ...


class AtomType(Enum): # FIX: write the actual values
    NONE:              AtomType # 0
    HbondAcceptor:     AtomType # 1
    HbondDonor:        AtomType # 2
    WeakHbondAcceptor: AtomType # 4
    WeakHbondDonor:    AtomType # 8
    PositiveCharge:    AtomType # 16
    NegativeCharge:    AtomType # 32
    CarbonylOxygen:    AtomType # 64
    CarbonylCarbon:    AtomType # 128
    Aromatic:          AtomType # 256
    Hydrophobic:       AtomType # 512
    XBondAcceptor:     AtomType # 1024
    XbondDonor:        AtomType # 2048
    IonicTypePartner:  AtomType # 4096
    DativeBondPartner: AtomType # 8192
    TransitionMetal:   AtomType # 16384
    IonicTypeMetal:    AtomType # 32768
    Invalid:           AtomType # 65536

    def __init__(self, value: int) -> None: ...

    def split(self) -> list[Self]: ...

    @property
    def value(self) -> str: ...

    def __index__ (self) -> int: ...
    def __int__   (self) -> int: ...
    def __hash__  (self) -> int: ...
    def __invert__(self) -> Self: ...

    def __and__(self, other: Self) -> Self: ...
    def __or__(self,  other: Self) -> Self: ...
    def __xor__(self, other: Self) -> Self: ...
    def __iand__(self, other: Self) -> Self: ...
    def __ior__(self,  other: Self) -> Self: ...
    def __ixor__(self, other: Self) -> Self: ...

    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...

    def __getstate__(self) -> int: ...
    def __setstate__(self, state: int) -> None: ...


class Flags:
    @staticmethod
    def has(types: AtomType, has_type: AtomType) -> bool: ...
    @staticmethod
    def has_any(types: AtomType, has_type: AtomType) -> bool: ...
    @staticmethod
    def all(types: AtomType, has_type: AtomType) -> bool: ...
    @staticmethod
    def any(types: AtomType, has_type: AtomType) -> bool: ...
    @staticmethod
    def none(types: AtomType, has_type: AtomType) -> bool: ...
    @staticmethod
    def empty(types: AtomType) -> bool: ...
    @staticmethod
    def split(types: AtomType) -> list[AtomType]: ...
    @staticmethod
    def get_enum_as_string(types: AtomType) -> str: ...
    @staticmethod
    def string_to_atom_type(s: str) -> AtomType: ...


def factorize_residues(resnames: Labels, resids: Indices, chainlabels: Labels) -> tuple[list[int], list[str], list[int], list[str]]: ...


class _EntityInterface(Protocol):
    def get_center(self) -> Point3D: ...
    def get_id(self) -> int: ...


_EI = TypeVar("_EI", bound=_EntityInterface)
class Entity(Generic[_EI]):
    def __init__(self, t: _EI) -> None: ...
    def get_center(self) -> Point3D: ...
    def get_id(self) -> int: ...

class Logger_:
    class LogLevel(Enum):
        Trace    = 0
        Debug    = 1
        Info     = 2
        Warn     = 3
        Error    = 4
        Critical = 5
        Off      = 6

    class FormatStyle(Enum):
        Simple   = 0
        Detailed = 1

    @staticmethod
    def get_instance() -> Logger_: ...

    def set_format(self,    style: FormatStyle) -> None: ...
    def set_log_level(self, level: LogLevel)    -> None: ...

    def log(self, level: LogLevel, message: str) -> None: ...



class PropertyKey(Enum):
    Names     = 0
    Indices   = 1
    Elements  = 2
    Positions = 3


class PropertyQueryLuni:
    def __init__(self) -> None: ...
    @overload
    def select(self, key: PropertyKey) -> PropertyQueryLuni: ...
    @overload
    def select(self, keys: list[PropertyKey]) -> PropertyQueryLuni: ...
    def properties(self) -> list[PropertyKey]: ...


class PropertyAnalyzerLuni:
    def __init__(self, query: PropertyQueryLuni) -> None: ...


class LuniResultType:
    def __init__(self) -> None: ...

    def has_property(self, key: PropertyKey) -> bool: ...
    @overload
    def get(self, key: Literal[PropertyKey.Names]) -> NDArray[np.str_]: ...
    @overload
    def get(self, key: Literal[PropertyKey.Indices]) -> NDArray[np.int32]: ...
    @overload
    def get(self, key: Literal[PropertyKey.Elements]) -> NDArray[np.str_]: ...
    @overload
    def get(self, key: Literal[PropertyKey.Positions]) -> NDArray[np.float64]: ...

    def get(self, key: PropertyKey) -> NDArray[np.str_] | NDArray[np.int32] | NDArray[np.float64]: ...

    @overload
    def __getitem__(self, key: Literal[PropertyKey.Names]) -> NDArray[np.str_]: ...
    @overload
    def __getitem__(self, key: Literal[PropertyKey.Indices]) -> NDArray[np.int32]: ...
    @overload
    def __getitem__(self, key: Literal[PropertyKey.Elements]) -> NDArray[np.str_]: ...
    @overload
    def __getitem__(self, key: Literal[PropertyKey.Positions]) -> NDArray[np.float64]: ...

    def __getitem__(self, key: PropertyKey) -> NDArray[np.str_] | NDArray[np.int32] | NDArray[np.float64]: ...

    def __contains__(self, key: PropertyKey) -> bool: ...


class LuniFileProcessor:
    def __init__(self, n_threads: int, analyzer: PropertyAnalyzerLuni, use_spinner: bool) -> None: ...
    def process_files(self, file_names: list[str]) -> None: ...
    def wait_for_completion(self) -> None: ...
    def get_result(self, file_name: str) -> LuniResultType | None: ...
    def get_all_results(self) -> dict[str, LuniResultType]: ...


@overload
def process_files(files: list[str], property_keys: list[PropertyKey], n_jobs: int, use_spinner: bool = True) -> dict[str, LuniResultType]: ...
@overload
def process_files(files: list[str], n_jobs: int, use_spinner: bool = True) -> dict[str, LahutaCPP]: ...


class AlignType(Enum):
    AA_3Di: int
    AA:     int
    _3Di:   int

class TMScoreThrMode(Enum):
    alignment: int
    query:     int
    target:    int

class SeqType(Enum):
    AminoAcid: int
    Nucleotide: int
    HMM: int

class FoldSeekOps:
    alignType: AlignType
    tmScoreThr: float
    tmScoreThrMode: TMScoreThrMode
    exactTMscore: bool
    lddtThr: float
    sortByStructureBits: bool
    alignmentMode: int
    alignmentOutputMode: int
    wrappedScoring: bool
    maxSeqLen: int
    compBiasCorrection: bool
    compBiasCorrectionScale: float
    scoreBias: float
    realign: bool
    correlationScoreWeight: float
    addBacktrace: int
    covThr: float
    covMode: int
    evalThr: int
    seqIdThr: float
    seqIdMode: bool
    alnLenThr: float
    chainNameMode: int
    maskBfactorThreshold: float
    altAlignment: int
    inputFormat: int
    gapOpen: int
    gapExtend: int

class PrefilterOptions:
    use_prefilter: bool
    alphabetSize: int
    maskMode: bool
    maskLowerCaseMode: bool
    maskProb: float
    kmerSize: int
    kmerThr: int
    spacedKmer: bool
    spacedKmerPattern: str
    takeOnlyBestKmer: bool
    querySeqType: SeqType
    targetSeqType: SeqType
    targetSearchMode: int
    sensitivity: float
    maxSeqLen: int
    diagonalScoring: int
    minDiagScoreThr: int
    aaBiasCorrection: bool
    aaBiasCorrectionScale: float
    covThr: float
    covMode: int
    maxResListLen: int

class ProcessingConfig:
    query_chunk_size: int
    target_chunk_size: int
    allow_self_ops: bool


FileList: TypeAlias = list[str]

class AlignerResults:
    query: SeqData
    target: SeqData
    results: list[Result]


class LahutaAlignerBase: ...
class LahutaAligner(LahutaAlignerBase):
    def __init__(self, ops: FoldSeekOps = ..., pf_ops: PrefilterOptions = ..., n_threads: int = 0) -> None: ...
    def run(self, query_files: FileList, target_files: FileList) -> None: ...
    def get_results(self) -> list[AlignerResults]: ...


class LahutaProcessor:
    def __init__(self, aligner: LahutaAlignerBase, config: ProcessingConfig) -> None: ...
    @overload
    def process(self, query_files: FileList, target_files: FileList) -> None: ...
    @overload
    def process(self, query_files: FileList) -> None: ...


class Result:
    dbKey: int
    score: int
    qcov: float
    dbcov: float
    seqId: float
    eval: float
    alnLength: int
    qStartPos: int
    qEndPos: int
    qLen: int
    dbStartPos: int
    dbEndPos: int
    dbLen: int
    queryOrfStartPos: int
    queryOrfEndPos: int
    dbOrfStartPos: int
    dbOrfEndPos: int
    backtrace: str

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self,
                 dbkey: int,
                 score: int,
                 qcov: float,
                 dbcov: float,
                 seqId: float,
                 eval: float,
                 alnLength: int,
                 qStartPos: int,
                 qEndPos: int,
                 qLen: int,
                 dbStartPos: int,
                 dbEndPos: int,
                 dbLen: int,
                 queryOrfStartPos: int,
                 queryOrfEndPos: int,
                 dbOrfStartPos: int,
                 dbOrfEndPos: int,
                 backtrace: str) -> None: ...

    @overload
    def __init__(self,
                 dbkey: int,
                 score: int,
                 qcov: float,
                 dbcov: float,
                 seqId: float,
                 eval: float,
                 alnLength: int,
                 qStartPos: int,
                 qEndPos: int,
                 qLen: int,
                 dbStartPos: int,
                 dbEndPos: int,
                 dbLen: int,
                 backtrace: str) -> None: ...

    def __init__(self, *args, **kwargs) -> None: ...

    @staticmethod
    def protein2nucl(backtrace: str, newBacktrace: str) -> None: ...


class Matcher:
    @staticmethod
    def compareHits(first: Result, second: Result) -> bool: ...

    @staticmethod
    def computeAlnLength(qStart: int, qEnd: int, dbStart: int, dbEnd: int) -> int: ...

    @staticmethod
    def compressAlignment(bt: str) -> str: ...

    @staticmethod
    def uncompressAlignment(cbt: str) -> str: ...

    @staticmethod
    def result_to_string(r: Result, compress_backtrace: bool = True) -> str: ...

    @staticmethod
    def results_to_string(vr: list[Result], compress_backtrace: bool = True) -> str: ...

    @staticmethod
    def resultToBuffer(result: Result, addBacktrace: bool, compress: bool = True, addOrfPosition: bool = False) -> str: ...


class SeqData:
    def __init__(self) -> None: ...
    def size(self) -> int: ...
    @property
    def x(self) -> NDArray[np.float64]: ...
    @property
    def y(self) -> NDArray[np.float64]: ...
    @property
    def z(self) -> NDArray[np.float64]: ...
    def __lt__(self, other: SeqData) -> bool: ...
    def map_3di(self, matrix: NDArray[np.float64], ops: FoldSeekOps) -> None: ...
    def map_aa(self, matrix: NDArray[np.float64], ops: FoldSeekOps) -> None: ...
    def build_sequence(self, matrix: NDArray[np.float64], ops: FoldSeekOps) -> None: ...

    Seq3Di: list[str]
    SeqAA: list[str]
    CaData: list[float]
    file_name: str
    chain_name: str


class TopologyMapper:
    class MappingType(IntEnum):
        Query = 1
        Target = 2

    def __init__(self, sd: SeqData, type: TopologyMapper.MappingType) -> None: ...  # noqa: A002
    def map(self, res: Result) -> None: ...
    def get_mapped_resid(self, atom_index: int) -> int | None: ...
    def get_entity_atoms(self, entity: EntityID) -> list[Atom]: ...
    def get_luni(self) -> LahutaCPP: ...


class ContactEquivKey:
    e1_mapped_ids: list[int]
    e2_mapped_ids: list[int]
    contact_type: InteractionType


class EquivalencyConfig:
    contact_resolution: bool
    contact_type: bool
    hbond_type: bool
    number_of_atoms: bool
    atom_name: bool
    element: bool
    resname: bool


class TopologicalEquivalency:
    def __init__(self, 
                 lm1: TopologyMapper,
                 lm2: TopologyMapper,
                 config: EquivalencyConfig | None = None) -> None: ...
    def evaluate(self, a1: Atom, a2: Atom) -> bool: ...
    def evaluate_atoms(self, a1: list[Atom], a2: list[Atom]) -> bool: ...
    def should_consider(self, atoms: list[Atom], mt: TopologyMapper.MappingType) -> bool: ...
    def is_mappable(self, mt: TopologyMapper.MappingType, atom: Atom) -> bool: ...
    def get_lm(self, type: TopologyMapper.MappingType) -> TopologyMapper: ...  # noqa: A002


class LahutaMapper:
    def __init__(self, query: SeqData, target: SeqData) -> None: ...
    def map(self, res: Result, config: EquivalencyConfig | None = None) -> None: ...
    def evaluate(self,
                 c1: Contacts,
                 m1: TopologyMapper.MappingType,
                 c2: Contacts,
                 m2: TopologyMapper.MappingType) -> int: ...
    def get_luni(self, type: TopologyMapper.MappingType) -> LahutaCPP: ...  # noqa: A002
