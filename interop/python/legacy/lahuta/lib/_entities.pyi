from collections.abc import Iterator
from enum import Enum
from typing import Callable, Generic, Protocol, TypeAlias, TypeVar, overload

from lahuta.lib._lahuta import Atom, AtomType, Flags, LahutaSystem, NSResults_, Point3D

# fmt: off

class FeatureGroup(Enum):
    NoneGroup       = 0
    QuaternaryAmine = 1
    TertiaryAmine   = 2
    Sulfonium       = 3
    SulfonicAcid    = 4
    Sulfate         = 5
    Phosphate       = 6
    Halocarbon      = 7
    Guanidine       = 8
    Acetamidine     = 9
    Carboxylate     = 10


class _EntityInterface(Protocol):
    @property
    def id(self) -> int: ...
    @property
    def center(self) -> Point3D: ...
    def has_atom(self, atom: Atom) -> bool: ...


class AtomEntity(_EntityInterface):
    type: AtomType
    atom: Atom

    @overload
    def __init__(self, luni: LahutaSystem, atom_index: int, atom_type: AtomType, center: Point3D, index: int) -> None: ...
    @overload
    def __init__(self, atom: Atom, atom_type: AtomType, center: Point3D, res_index: int) -> None: ...

    def get_data(self) -> Atom: ...


class RingEntity(_EntityInterface):
    atoms: list[Atom]

    @overload
    def __init__(self, luni: LahutaSystem, atom_ids: list[int], center: Point3D, normal: Point3D, index: int) -> None: ...
    @overload
    def __init__(self, atoms: list[Atom], center: Point3D, normal: Point3D, index: int) -> None: ...

    def get_data(self) -> list[Atom]: ...
    def get_atom_ids(self) -> list[int]: ...

    @property
    def normal(self) -> Point3D: ...


class GroupEntity(_EntityInterface):
    type: AtomType
    group: FeatureGroup
    atoms: list[Atom]

    @overload
    def __init__(self, luni: LahutaSystem, atom_ids: list[int],  atom_type: AtomType, group: FeatureGroup, center: Point3D, index: int) -> None: ...
    @overload
    def __init__(self, aatoms: list[Atom], tom_type: AtomType, group: FeatureGroup, center: Point3D, index: int) -> None: ...

    def get_data(self) -> list[Atom]: ...
    def get_atom_ids(self) -> list[int]: ...


_Entity = TypeVar("_Entity", AtomEntity, RingEntity, GroupEntity)
class EntityCollection(Generic[_Entity]):
    data: list[_Entity]

    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, entities: list[_Entity]) -> None: ...

    def add_entity(self, entity: _Entity) -> None: ...
    def get_entities(self) -> list[_Entity]: ...
    def get_positions(self) -> list[Point3D]: ...
    def get_atom_ids(self) -> list[int]: ...
    def get_atoms(self) -> list[Atom]: ...
    def get_size(self) -> int: ...
    def __getitem__(self, index: int) -> _Entity: ...
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[_Entity]: ...


AtomTypeCheckFunc: TypeAlias = Callable[[AtomType, AtomType], bool]

class AtomEntityCollection(EntityCollection[AtomEntity]):
    @staticmethod
    def filter(luni: LahutaSystem, type_: AtomType, check_func: AtomTypeCheckFunc = Flags.has_any) -> AtomEntityCollection: ...

class RingEntityCollection(EntityCollection[RingEntity]):
    rings: list[RingEntity]

class GroupEntityCollection(EntityCollection[GroupEntity]):
    @staticmethod
    def filter(luni: LahutaSystem, type_: AtomType, check_func: AtomTypeCheckFunc = Flags.has_any) -> GroupEntityCollection: ...


# Works, but the type inference is not as good as with exhaustive (and painfully annoying) overloads.
# _EC1 = TypeVar("_EC1", AtomEntityCollection, RingEntityCollection, GroupEntityCollection)
# _EC2 = TypeVar("_EC2", AtomEntityCollection, RingEntityCollection, GroupEntityCollection)
# class EntityNeighborSearch:
#     @staticmethod
#     def search(a: _EC1, b: _EC2, cutoff: float) -> NSResults_:
#         ...

class EntityNeighborSearch:
    @overload
    @staticmethod
    def search(a: AtomEntityCollection,  b: AtomEntityCollection,  cutoff: float) -> NSResults_: ...

    @overload
    @staticmethod
    def search(a: AtomEntityCollection,  b: RingEntityCollection,  cutoff: float) -> NSResults_: ...

    @overload
    @staticmethod
    def search(a: AtomEntityCollection,  b: GroupEntityCollection, cutoff: float) -> NSResults_: ...

    @overload
    @staticmethod
    def search(a: RingEntityCollection,  b: AtomEntityCollection,  cutoff: float) -> NSResults_: ...

    @overload
    @staticmethod
    def search(a: RingEntityCollection,  b: RingEntityCollection,  cutoff: float) -> NSResults_: ...

    @overload
    @staticmethod
    def search(a: RingEntityCollection,  b: GroupEntityCollection, cutoff: float) -> NSResults_: ...

    @overload
    @staticmethod
    def search(a: GroupEntityCollection, b: AtomEntityCollection,  cutoff: float) -> NSResults_: ...

    @overload
    @staticmethod
    def search(a: GroupEntityCollection, b: RingEntityCollection,  cutoff: float) -> NSResults_: ...

    @overload
    @staticmethod
    def search(a: GroupEntityCollection, b: GroupEntityCollection, cutoff: float) -> NSResults_: ...
