"""
Bindings for RDKit core types
"""
from __future__ import annotations
import numpy
import typing
__all__: list[str] = ['AROMATIC', 'Atom', 'AtomMonomerInfo', 'AtomPDBResidueInfo', 'BEGINDASH', 'BEGINWEDGE', 'Bond', 'BondDir', 'BondStereo', 'BondType', 'Conformer', 'DOUBLE', 'EITHERDOUBLE', 'ENDDOWNRIGHT', 'ENDUPRIGHT', 'NONE', 'Point3D', 'RWMol', 'SINGLE', 'STEREOANY', 'STEREOATROPCCW', 'STEREOATROPCW', 'STEREOCIS', 'STEREOE', 'STEREONONE', 'STEREOTRANS', 'STEREOZ', 'TRIPLE', 'UNKNOWN', 'UNSPECIFIED', 'hasNonZeroZCoords', 'pAtomPDBResidueInfo']
class Atom:
    class HybridizationType:
        OTHER: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.OTHER: 8>
        S: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.S: 1>
        SP: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.SP: 2>
        SP2: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.SP2: 3>
        SP2D: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.SP2D: 5>
        SP3: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.SP3: 4>
        SP3D: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.SP3D: 6>
        SP3D2: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.SP3D2: 7>
        UNSPECIFIED: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.UNSPECIFIED: 0>
        __members__: typing.ClassVar[dict[str, Atom.HybridizationType]]  # value = {'UNSPECIFIED': <HybridizationType.UNSPECIFIED: 0>, 'S': <HybridizationType.S: 1>, 'SP': <HybridizationType.SP: 2>, 'SP2': <HybridizationType.SP2: 3>, 'SP3': <HybridizationType.SP3: 4>, 'SP2D': <HybridizationType.SP2D: 5>, 'SP3D': <HybridizationType.SP3D: 6>, 'SP3D2': <HybridizationType.SP3D2: 7>, 'OTHER': <HybridizationType.OTHER: 8>}
        def __eq__(self, other: typing.Any) -> bool:
            ...
        def __getstate__(self) -> int:
            ...
        def __hash__(self) -> int:
            ...
        def __index__(self) -> int:
            ...
        def __init__(self, value: int) -> None:
            ...
        def __int__(self) -> int:
            ...
        def __ne__(self, other: typing.Any) -> bool:
            ...
        def __repr__(self) -> str:
            ...
        def __setstate__(self, state: int) -> None:
            ...
        def __str__(self) -> str:
            ...
        @property
        def name(self) -> str:
            ...
        @property
        def value(self) -> int:
            ...
    OTHER: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.OTHER: 8>
    S: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.S: 1>
    SP: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.SP: 2>
    SP2: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.SP2: 3>
    SP2D: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.SP2D: 5>
    SP3: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.SP3: 4>
    SP3D: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.SP3D: 6>
    SP3D2: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.SP3D2: 7>
    UNSPECIFIED: typing.ClassVar[Atom.HybridizationType]  # value = <HybridizationType.UNSPECIFIED: 0>
    def calcExplicitValence(self, strict: bool = True) -> int:
        ...
    def calcImplicitValence(self, strict: bool = True) -> int:
        ...
    def getAtomicNum(self) -> int:
        ...
    def getDegree(self) -> int:
        ...
    def getExplicitValence(self) -> int:
        ...
    def getFormalCharge(self) -> int:
        ...
    def getHybridization(self) -> ...:
        ...
    def getIdx(self) -> int:
        ...
    def getImplicitValence(self) -> int:
        ...
    def getIsAromatic(self) -> bool:
        ...
    def getMass(self) -> float:
        ...
    def getMonomerInfo(self) -> AtomPDBResidueInfo:
        """
        PDB residue annotation for this atom, or None if absent or non-PDB type.
        """
    def getNumExplicitHs(self) -> int:
        ...
    def getNumImplicitHs(self) -> int:
        ...
    def getSymbol(self) -> str:
        ...
    def getTotalNumHs(self, include_neighbors: bool = False) -> int:
        ...
    def getTotalValence(self) -> int:
        ...
    def setAtomicNum(self, arg0: int) -> None:
        ...
class AtomMonomerInfo:
    """
    Per-atom monomer metadata; owned by the Atom.
    """
    class AtomMonomerType:
        OTHER: typing.ClassVar[AtomMonomerInfo.AtomMonomerType]  # value = <AtomMonomerType.OTHER: 2>
        PDBRESIDUE: typing.ClassVar[AtomMonomerInfo.AtomMonomerType]  # value = <AtomMonomerType.PDBRESIDUE: 1>
        UNKNOWN: typing.ClassVar[AtomMonomerInfo.AtomMonomerType]  # value = <AtomMonomerType.UNKNOWN: 0>
        __members__: typing.ClassVar[dict[str, AtomMonomerInfo.AtomMonomerType]]  # value = {'UNKNOWN': <AtomMonomerType.UNKNOWN: 0>, 'PDBRESIDUE': <AtomMonomerType.PDBRESIDUE: 1>, 'OTHER': <AtomMonomerType.OTHER: 2>}
        def __eq__(self, other: typing.Any) -> bool:
            ...
        def __getstate__(self) -> int:
            ...
        def __hash__(self) -> int:
            ...
        def __index__(self) -> int:
            ...
        def __init__(self, value: int) -> None:
            ...
        def __int__(self) -> int:
            ...
        def __ne__(self, other: typing.Any) -> bool:
            ...
        def __repr__(self) -> str:
            ...
        def __setstate__(self, state: int) -> None:
            ...
        def __str__(self) -> str:
            ...
        @property
        def name(self) -> str:
            ...
        @property
        def value(self) -> int:
            ...
    OTHER: typing.ClassVar[AtomMonomerInfo.AtomMonomerType]  # value = <AtomMonomerType.OTHER: 2>
    PDBRESIDUE: typing.ClassVar[AtomMonomerInfo.AtomMonomerType]  # value = <AtomMonomerType.PDBRESIDUE: 1>
    UNKNOWN: typing.ClassVar[AtomMonomerInfo.AtomMonomerType]  # value = <AtomMonomerType.UNKNOWN: 0>
    def copy(self) -> AtomMonomerInfo:
        """
        Deep-copy the monomer info
        """
    def getMonomerType(self) -> AtomMonomerInfo.AtomMonomerType:
        """
        Kind of monomer annotation
        """
    def getName(self) -> str:
        """
        Monomer name (e.g., atom label)
        """
    def setMonomerType(self, arg0: AtomMonomerInfo.AtomMonomerType) -> None:
        ...
    def setName(self, name: str) -> None:
        """
        Set the monomer name
        """
class AtomPDBResidueInfo(AtomMonomerInfo):
    """
    PDB residue-level annotation for an atom (serial, resname, chain, etc.).
    """
    def __init__(self, atomName: str, serialNumber: int, residueName: str, residueNumber: int, chainId: str, altLoc: str = '', insertionCode: str = '', occupancy: float = 1.0, tempFactor: float = 0.0, isHet: bool = False, secondaryStructure: int = 0, segmentNumber: int = 0) -> None:
        """
        Create a PDB residue annotation (requires 5 fields: atomName, serialNumber, residueName, residueNumber, chainId).
        """
    def getAltLoc(self) -> str:
        """
        Get the alternate location identifier (altLoc) for the atom
        """
    def getChainId(self) -> str:
        """
        Get the chain identifier for the atom
        """
    def getInsertionCode(self) -> str:
        """
        Get the insertion code for the atom
        """
    def getIsHeteroAtom(self) -> bool:
        """
        Check if the atom is a hetero atom (non-standard residue)
        """
    def getOccupancy(self) -> float:
        """
        Get the occupancy value for the atom
        """
    def getResidueIndex(self) -> int:
        """
        Index of the residue in source (non-PDB standard)
        """
    def getResidueName(self) -> str:
        """
        Get the residue name
        """
    def getResidueNumber(self) -> int:
        """
        Get the residue number for the atom
        """
    def getSecondaryStructure(self) -> int:
        """
        Get the secondary structure identifier for the atom
        """
    def getSegmentNumber(self) -> int:
        """
        Get the segment number for the atom
        """
    def getSerialNumber(self) -> int:
        """
        Get the PDB serial number of the atom
        """
    def getTempFactor(self) -> float:
        """
        Get the temperature factor (B-factor) for the atom
        """
    def setAltLoc(self, altLoc: str) -> None:
        """
        Set the alternate location identifier (altLoc) for the atom
        """
    def setChainId(self, chainId: str) -> None:
        """
        Set the chain identifier for the atom
        """
    def setInsertionCode(self, insertionCode: str) -> None:
        """
        Set the insertion code for the atom
        """
    def setIsHeteroAtom(self, isHeteroAtom: bool) -> None:
        """
        Set whether the atom is a hetero atom (non-standard residue)
        """
    def setOccupancy(self, occupancy: float) -> None:
        """
        Set the occupancy value for the atom
        """
    def setResidueIndex(self, residueIndex: int) -> None:
        """
        Set the residue index for the atom (non-PDB standard)
        """
    def setResidueName(self, residueName: str) -> None:
        """
        Set the residue name for the atom
        """
    def setResidueNumber(self, residueNumber: int) -> None:
        """
        Set the residue number for the atom
        """
    def setSecondaryStructure(self, secondaryStructure: int) -> None:
        """
        Set the secondary structure identifier for the atom
        """
    def setSegmentNumber(self, segmentNumber: int) -> None:
        """
        Set the segment number for the atom
        """
    def setSerialNumber(self, serialNumber: int) -> None:
        """
        Set the PDB serial number of the atom
        """
    def setTempFactor(self, tempFactor: float) -> None:
        """
        Set the temperature factor (B-factor) for the atom
        """
class Bond:
    """
    RDKit bond object (owned by a molecule)
    """
    def __str__(self) -> str:
        ...
    def getBeginAtom(self) -> Atom:
        """
        Begin atom (by reference)
        """
    def getBeginAtomIdx(self) -> int:
        ...
    def getBondDir(self) -> BondDir:
        """
        Directional annotation
        """
    def getBondType(self) -> BondType:
        """
        Bond type (order)
        """
    def getBondTypeAsDouble(self) -> float:
        """
        Bond order as float (e.g., 1.0, 1.5, 2.0)
        """
    def getEndAtom(self) -> Atom:
        """
        End atom (by reference)
        """
    def getEndAtomIdx(self) -> int:
        ...
    def getIdx(self) -> int:
        """
        Index within owning molecule
        """
    def getIsAromatic(self) -> bool:
        ...
    def getIsConjugated(self) -> bool:
        ...
    def getOtherAtomIdx(self, this_idx: int) -> int:
        ...
    def getStereo(self) -> BondStereo:
        """
        Stereo annotation
        """
    def getStereoAtoms(self) -> list:
        """
        Indices of stereo reference atoms (size 0 or 2)
        """
    def hasOwningMol(self) -> bool:
        ...
    def setBondDir(self, arg0: BondDir) -> None:
        ...
    def setBondType(self, arg0: BondType) -> None:
        ...
    def setIsAromatic(self, arg0: bool) -> None:
        ...
    def setIsConjugated(self, arg0: bool) -> None:
        ...
    def setStereo(self, arg0: BondStereo) -> None:
        ...
    def setStereoAtoms(self, begin_neighbor_idx: int, end_neighbor_idx: int) -> None:
        """
        Set neighboring atoms used as reference for cis/trans
        """
class BondDir:
    BEGINDASH: typing.ClassVar[BondDir]  # value = <BondDir.BEGINDASH: 2>
    BEGINWEDGE: typing.ClassVar[BondDir]  # value = <BondDir.BEGINWEDGE: 1>
    EITHERDOUBLE: typing.ClassVar[BondDir]  # value = <BondDir.EITHERDOUBLE: 5>
    ENDDOWNRIGHT: typing.ClassVar[BondDir]  # value = <BondDir.ENDDOWNRIGHT: 3>
    ENDUPRIGHT: typing.ClassVar[BondDir]  # value = <BondDir.ENDUPRIGHT: 4>
    NONE: typing.ClassVar[BondDir]  # value = <BondDir.NONE: 0>
    UNKNOWN: typing.ClassVar[BondDir]  # value = <BondDir.UNKNOWN: 6>
    __members__: typing.ClassVar[dict[str, BondDir]]  # value = {'NONE': <BondDir.NONE: 0>, 'BEGINWEDGE': <BondDir.BEGINWEDGE: 1>, 'BEGINDASH': <BondDir.BEGINDASH: 2>, 'ENDDOWNRIGHT': <BondDir.ENDDOWNRIGHT: 3>, 'ENDUPRIGHT': <BondDir.ENDUPRIGHT: 4>, 'EITHERDOUBLE': <BondDir.EITHERDOUBLE: 5>, 'UNKNOWN': <BondDir.UNKNOWN: 6>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class BondStereo:
    STEREOANY: typing.ClassVar[BondStereo]  # value = <BondStereo.STEREOANY: 1>
    STEREOATROPCCW: typing.ClassVar[BondStereo]  # value = <BondStereo.STEREOATROPCCW: 7>
    STEREOATROPCW: typing.ClassVar[BondStereo]  # value = <BondStereo.STEREOATROPCW: 6>
    STEREOCIS: typing.ClassVar[BondStereo]  # value = <BondStereo.STEREOCIS: 4>
    STEREOE: typing.ClassVar[BondStereo]  # value = <BondStereo.STEREOE: 3>
    STEREONONE: typing.ClassVar[BondStereo]  # value = <BondStereo.STEREONONE: 0>
    STEREOTRANS: typing.ClassVar[BondStereo]  # value = <BondStereo.STEREOTRANS: 5>
    STEREOZ: typing.ClassVar[BondStereo]  # value = <BondStereo.STEREOZ: 2>
    __members__: typing.ClassVar[dict[str, BondStereo]]  # value = {'STEREONONE': <BondStereo.STEREONONE: 0>, 'STEREOANY': <BondStereo.STEREOANY: 1>, 'STEREOZ': <BondStereo.STEREOZ: 2>, 'STEREOE': <BondStereo.STEREOE: 3>, 'STEREOCIS': <BondStereo.STEREOCIS: 4>, 'STEREOTRANS': <BondStereo.STEREOTRANS: 5>, 'STEREOATROPCW': <BondStereo.STEREOATROPCW: 6>, 'STEREOATROPCCW': <BondStereo.STEREOATROPCCW: 7>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class BondType:
    AROMATIC: typing.ClassVar[BondType]  # value = <BondType.AROMATIC: 12>
    DOUBLE: typing.ClassVar[BondType]  # value = <BondType.DOUBLE: 2>
    SINGLE: typing.ClassVar[BondType]  # value = <BondType.SINGLE: 1>
    TRIPLE: typing.ClassVar[BondType]  # value = <BondType.TRIPLE: 3>
    UNSPECIFIED: typing.ClassVar[BondType]  # value = <BondType.UNSPECIFIED: 0>
    __members__: typing.ClassVar[dict[str, BondType]]  # value = {'UNSPECIFIED': <BondType.UNSPECIFIED: 0>, 'SINGLE': <BondType.SINGLE: 1>, 'DOUBLE': <BondType.DOUBLE: 2>, 'TRIPLE': <BondType.TRIPLE: 3>, 'AROMATIC': <BondType.AROMATIC: 12>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class Conformer:
    @typing.overload
    def __init__(self) -> None:
        """
        Create an empty conformer
        """
    @typing.overload
    def __init__(self, numAtoms: int) -> None:
        """
        Create with N atoms initialized to (0,0,0)
        """
    def getAtomPos(self, idx: int) -> ...:
        """
        Get atom position by index
        """
    def getId(self) -> int:
        """
        Conformer identifier
        """
    def getNumAtoms(self) -> int:
        """
        Number of atom positions
        """
    def getPositions(self) -> numpy.ndarray[typing.Any, numpy.dtype[typing.Any]]:
        """
        Atomic positions as a zero-copy NumPy view
        """
    def hasOwningMol(self) -> bool:
        """
        Whether this conformer is attached to a molecule
        """
    def is3D(self) -> bool:
        """
        Whether positions are 3D
        """
    def reserve(self, size: int) -> None:
        """
        Reserve capacity for positions
        """
    def resize(self, size: int) -> None:
        """
        Resize internal position vector
        """
    def set3D(self, arg0: bool) -> None:
        ...
    def setAllAtomPositions(self, positions: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]) -> None:
        """
        Replace all positions from numpy array (N,3)
        """
    @typing.overload
    def setAtomPos(self, idx: int, pos: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]) -> None:
        """
        Set atom position from numpy array of shape (3,)
        """
    @typing.overload
    def setAtomPos(self, idx: int, pos: ...) -> None:
        """
        Set atom position from Point3D
        """
    def setId(self, arg0: int) -> None:
        ...
class Point3D:
    x: float
    y: float
    z: float
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: float, arg1: float, arg2: float) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]) -> None:
        ...
    def __str__(self) -> str:
        ...
class RWMol:
    def __getitem__(self, idx: int) -> Atom:
        ...
    def __init__(self) -> None:
        ...
    def __iter__(self) -> typing.Any:
        ...
    @typing.overload
    def addAtom(self, update_label: bool = True) -> int:
        """
        Add a new (empty) atom; returns new atom index
        """
    @typing.overload
    def addAtom(self, atomic_number: int, update_label: bool = True) -> int:
        """
        Add a new atom with a given atomic number; returns new atom index
        """
    @typing.overload
    def addBond(self, begin_idx: int, end_idx: int) -> int:
        """
        Add a bond with unspecified order. Returns new number of bonds
        """
    @typing.overload
    def addBond(self, begin_idx: int, end_idx: int, order: BondType) -> int:
        """
        Add a bond; returns new number of bonds
        """
    def addConformer(self, conformer: Conformer, assign_id: bool = False) -> int:
        """
        Add a conformer (copied) to the molecule; returns its id
        """
    def atomBondObjects(self, idx: int) -> list:
        """
        Bond objects incident to the given atom (by reference)
        """
    def atomBonds(self, idx: int) -> list:
        """
        Bond endpoint index pairs for bonds incident to the given atom
        """
    def atomNeighbors(self, idx: int) -> list:
        """
        Neighbor atom indices of the given atom
        """
    def atoms(self) -> list:
        """
        List of Atom objects (by reference)
        """
    def bondPairs(self) -> list:
        """
        List of (begin_idx, end_idx) for all bonds
        """
    def bonds(self) -> list:
        """
        List of Bond objects (by reference)
        """
    def clearComputedProps(self, include_rings: bool = True) -> None:
        """
        Clear computed properties and optionally ring info
        """
    def clearConformers(self) -> None:
        """
        Remove all conformers
        """
    def debugMol(self) -> str:
        """
        Return a debug string for the molecule
        """
    def getAtomWithIdx(self, idx: int) -> Atom:
        """
        Get atom by index (by reference)
        """
    def getBondBetweenAtoms(self, i: int, j: int) -> Bond:
        """
        Get bond between atom indices i and j, or None
        """
    def getBondWithIdx(self, idx: int) -> Bond:
        """
        Get bond by index (by reference)
        """
    def getConformer(self, id: int = -1) -> Conformer:
        """
        Get a conformer by id (or the first if id < 0)
        """
    def getConnectedComponents(self) -> tuple:
        """
        Return (num_components, component_ids)
        """
    @typing.overload
    def getNumAtoms(self) -> int:
        """
        Total number of atoms (explicit only)
        """
    @typing.overload
    def getNumAtoms(self, only_explicit: bool) -> int:
        """
        Number of atoms; include implicit Hs when only_explicit is False
        """
    @typing.overload
    def getNumBonds(self) -> int:
        """
        Number of bonds (heavy atoms only)
        """
    @typing.overload
    def getNumBonds(self, only_heavy: bool = True) -> int:
        """
        Number of bonds; include H-bonds when only_heavy is False
        """
    def getNumConformers(self) -> int:
        """
        Number of conformers in the molecule
        """
    def getPositions(self, conf_id: int = -1) -> numpy.ndarray[typing.Any, numpy.dtype[typing.Any]]:
        """
        Get atomic positions as a zero-copy NumPy view for the specified conformer
        """
    def removeAtom(self, idx: int) -> None:
        """
        Remove an atom by index
        """
    def removeBond(self, begin_idx: int, end_idx: int) -> None:
        """
        Remove a bond between two atom indices
        """
    def removeConformer(self, id: int) -> None:
        """
        Remove conformer with a given id
        """
    def updatePropertyCache(self, strict: bool = True) -> None:
        """
        Recompute atom/bond property caches
        """
class pAtomPDBResidueInfo(AtomPDBResidueInfo):
    """
    Pooled AtomPDBResidueInfo optimized for Lahuta's object pool.
    """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, atomName: str, serialNumber: int, residueName: str, residueNumber: int) -> None:
        """
        Construct with common PDB fields. Other fields default.
        """
    def initialize(self, atomName: str, serialNumber: int, residueName: str, residueNumber: int) -> None:
        """
        Reinitialize fields without allocating a new object.
        """
    def resetState(self) -> None:
        """
        Reset fields to defaults (pool-friendly).
        """
def hasNonZeroZCoords(conformer: Conformer) -> bool:
    """
    Return True if any z coordinate has non-zero magnitude.
    """
AROMATIC: BondType  # value = <BondType.AROMATIC: 12>
BEGINDASH: BondDir  # value = <BondDir.BEGINDASH: 2>
BEGINWEDGE: BondDir  # value = <BondDir.BEGINWEDGE: 1>
DOUBLE: BondType  # value = <BondType.DOUBLE: 2>
EITHERDOUBLE: BondDir  # value = <BondDir.EITHERDOUBLE: 5>
ENDDOWNRIGHT: BondDir  # value = <BondDir.ENDDOWNRIGHT: 3>
ENDUPRIGHT: BondDir  # value = <BondDir.ENDUPRIGHT: 4>
NONE: BondDir  # value = <BondDir.NONE: 0>
SINGLE: BondType  # value = <BondType.SINGLE: 1>
STEREOANY: BondStereo  # value = <BondStereo.STEREOANY: 1>
STEREOATROPCCW: BondStereo  # value = <BondStereo.STEREOATROPCCW: 7>
STEREOATROPCW: BondStereo  # value = <BondStereo.STEREOATROPCW: 6>
STEREOCIS: BondStereo  # value = <BondStereo.STEREOCIS: 4>
STEREOE: BondStereo  # value = <BondStereo.STEREOE: 3>
STEREONONE: BondStereo  # value = <BondStereo.STEREONONE: 0>
STEREOTRANS: BondStereo  # value = <BondStereo.STEREOTRANS: 5>
STEREOZ: BondStereo  # value = <BondStereo.STEREOZ: 2>
TRIPLE: BondType  # value = <BondType.TRIPLE: 3>
UNKNOWN: BondDir  # value = <BondDir.UNKNOWN: 6>
UNSPECIFIED: BondType  # value = <BondType.UNSPECIFIED: 0>
