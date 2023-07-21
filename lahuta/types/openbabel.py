from typing import Any, List, Protocol


class ObSmartPatternType(Protocol):
    def Init(self, smarts: str) -> None:
        ...

    def Match(self, mol: Any) -> None:
        ...

    def GetMapList(self) -> List[Any]:
        ...


# pylint: disable=C0103
class OBSmartsPatternWrapper:
    def __init__(self, ob_smarts_pattern: ObSmartPatternType):
        self.ob_smarts_pattern = ob_smarts_pattern

    def Init(self, smarts: str) -> None:
        return self.ob_smarts_pattern.Init(smarts)

    def Match(self, mol: Any) -> None:
        return self.ob_smarts_pattern.Match(mol)

    def GetMapList(self) -> List[Any]:
        return self.ob_smarts_pattern.GetMapList()


# class


# pylint: disable=C0103
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
