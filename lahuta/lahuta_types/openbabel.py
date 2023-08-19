"""Provides a typed interface and typed wrappers for OpenBabel, a chemical toolbox
for chemical data. This module allows for static typing support when working with OpenBabel objects.

The module is comprised of several classes that correspond to various aspects of OpenBabel's functionality. 
Each class represents a specific object type or command in OpenBabel, with typed method interfaces.

The goal of this module is to allow for more rigorous static type checking and error detection when using OpenBabel.
However, in the future, OpenBabel may be replaced with RDKit, a more widely used, actively maintained, 
and more community-engaged library. If this transition occurs, hopefully RDKit will have native typing support,
eliminating the need for this module.

Classes:
    ObSmartPatternType(Protocol), OBSmartsPatternWrapper: Types and wrapper for OpenBabel's smarts pattern.
    MolType(Protocol), MolTypeWrapper: Types and wrapper for OpenBabel's molecule object.
    MolAtomType(Protocol): Types for OpenBabel's atom within a molecule.
    MolResType(Protocol): Types for OpenBabel's residue within a molecule.
    MolBond(Protocol): Types for OpenBabel's bond within a molecule.
    BondIterator(Protocol), BondIterable(Protocol): Iterator type for OpenBabel's bond object.
    ObRingType(Protocol): Types for OpenBabel's ring within a molecule.
    ObVector3Type(Protocol), ObVector3Wrapper: Types and wrapper for OpenBabel's 3D vector object.

Example:
    # Assuming appropriate OpenBabel objects are available
    ob_pattern = OBSmartsPatternWrapper(ob_smarts_pattern)
    ob_pattern.Init('c1ccccc1')
    ob_pattern.Match(mol)
    ob_pattern.GetMapList()

"""

from typing import Any, Iterator, List, Protocol

# pylint: disable=C0116, C0115, C0103


class ObSmartPatternType(Protocol):
    def Init(self, smarts: str) -> None:
        ...

    def Match(self, mol: Any) -> None:
        ...

    def GetMapList(self) -> List[Any]:
        ...


class OBSmartsPatternWrapper:
    def __init__(self, ob_smarts_pattern: ObSmartPatternType):
        self.ob_smarts_pattern = ob_smarts_pattern

    def Init(self, smarts: str) -> None:
        """Initialize the SMARTS pattern."""
        return self.ob_smarts_pattern.Init(smarts)

    def Match(self, mol: Any) -> None:
        """Match the SMARTS pattern to a molecule."""
        return self.ob_smarts_pattern.Match(mol)

    def GetMapList(self) -> List[Any]:
        """Get the map list of the SMARTS pattern."""
        return self.ob_smarts_pattern.GetMapList()


# class


class MolType(Protocol):
    def NumAtoms(self) -> int:
        ...

    def GetAtom(self, index: int) -> Any:
        ...

    # def GetId(self) -> Any:
    #     ...

    def GetAtomById(self, index: int) -> Any:
        ...

    def NewAtom(self, index: int) -> "MolAtomType":
        ...

    def NewBond(self) -> Any:
        ...

    def NewResidue(self) -> "MolResType":
        ...

    def ConnectTheDots(self) -> None:
        ...

    def PerceiveBondOrders(self) -> None:
        ...

    def BeginModify(self) -> Any:
        ...

    def EndModify(self, flag: bool = True) -> None:
        ...

    def SetAromaticPerceived(self) -> None:
        ...

    def SetAtomTypesPerceived(self) -> None:
        ...

    def SetChiralityPerceived(self) -> None:
        ...

    def SetRingTypesPerceived(self) -> None:
        ...

    def SetPartialChargesPerceived(self) -> None:
        ...

    def SetChainsPerceived(self) -> None:
        ...

    def NumBonds(self) -> int:
        ...

    def GetSSSR(self) -> Iterator["ObRingType"]:
        ...


class MolTypeWrapper:
    def __init__(self, mol: MolType):
        self.mol = mol

    def GetId(self) -> int:
        ...


class MolAtomType(Protocol):
    def GetId(self) -> int:
        ...

    def GetIdx(self) -> int:
        ...

    def GetAtomicNum(self) -> int:
        ...

    def SetType(self, atom_name: str) -> None:
        ...

    def SetPartialCharge(self, charge: float) -> None:
        ...

    def SetVector(self, x: float, y: float, z: float) -> None:
        ...

    def SetFormalCharge(self, charge: int) -> None:
        ...

    def SetAtomicNum(self, atomic_num: int) -> None:
        ...


class MolResType(Protocol):
    def AddAtom(self, atom: MolAtomType) -> None:
        ...

    def SetHetAtom(self, atom: MolAtomType, is_het: bool) -> None:
        ...

    def SetSerialNum(self, atom: MolAtomType, serial_num: int) -> None:
        ...

    def GetSerialNum(self, atom: MolAtomType) -> int:
        ...

    def SetChainNum(self, chain_num: int) -> None:
        ...

    def SetNum(self, num: str) -> None:
        ...

    def SetName(self, name: str) -> None:
        ...


class MolBond(Protocol):
    def GetBeginAtomIdx(self) -> int:
        ...

    def GetEndAtomIdx(self) -> int:
        ...


class BondIterator(Protocol):
    def __next__(self) -> MolBond:
        ...

    def __iter__(self) -> "BondIterator":
        ...


class BondIterable(Protocol):
    def __iter__(self) -> BondIterator:
        ...


class ObRingType(Protocol):
    def IsAromatic(self) -> bool:
        ...

    def findCenterAndNormal(self, center: "ObVector3Type", normal: "ObVector3Type", point: "ObVector3Type") -> None:
        ...

    @property
    def _path(self) -> List[int]:
        ...


class ObVector3Type(Protocol):
    def GetX(self) -> float:
        ...

    def GetY(self) -> float:
        ...

    def GetZ(self) -> float:
        ...


class ObVector3Wrapper:
    def __init__(self, center: ObVector3Type, normal: ObVector3Type) -> None:
        self.center = center
        self.normal = normal
