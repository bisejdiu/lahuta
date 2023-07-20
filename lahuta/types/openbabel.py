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


# pylint: disable=C0103
class MolType(Protocol):
    def NumAtoms(self) -> int:
        ...

    def GetAtom(self, index: int) -> Any:
        ...

    def GetId(self) -> int:
        ...
