from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, Iterator, List, Optional, Union

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from lahuta.core.loaders import GemmiLoader, TopologyLoader

from lahuta.types.mdanalysis import AtomGroupType, UniverseType


class Atoms:
    dtype = np.dtype(
        {
            "names": ["name", "id", "element", "type"],
            "formats": ["<U10", "int", "<U10", "<U10"],
        }
    )

    def __init__(self) -> None:
        self._data: NDArray[Any] = np.empty(0, dtype=self.dtype)
        self._coordinates = np.zeros((0, 3), dtype=np.float_)

    @classmethod
    def from_gemmi(cls, gemmi_block: Dict[str, Any]) -> "Atoms":
        cls_instance = cls.__new__(cls)

        # Create structured array
        label_atom_id: List[str] = gemmi_block["label_atom_id"]
        data = np.empty(len(label_atom_id), dtype=cls_instance.dtype)
        data["name"] = np.array(label_atom_id)
        data["id"] = np.array(gemmi_block.get("id"), dtype=int)
        data["element"] = np.array(gemmi_block.get("type_symbol"))
        data["type"] = np.array(gemmi_block.get("type_symbol"))

        cls_instance._data = data
        cls_instance._coordinates = np.zeros((0, 3))

        return cls_instance

    @classmethod
    def from_mda(cls, uv: UniverseType) -> "Atoms":
        cls_instance = cls.__new__(cls)
        data: NDArray[Any] = np.empty(len(uv.atoms), dtype=cls_instance.dtype)
        data["name"] = uv.atoms.names
        data["id"] = uv.atoms.ids
        data["element"] = uv.atoms.elements
        data["type"] = uv.atoms.types

        cls_instance._data = data
        cls_instance._coordinates = uv.atoms.positions

        return cls_instance

    @property
    def names(self) -> NDArray[np.str_]:
        return self._data["name"]

    @property
    def types(self) -> NDArray[np.str_]:
        return self._data["type"]

    @property
    def ids(self) -> NDArray[np.int32]:
        return self._data["id"]

    @property
    def elements(self) -> NDArray[np.str_]:
        return self._data["element"]

    @property
    def coordinates(self) -> NDArray[np.float_]:
        return self._coordinates

    @coordinates.setter
    def coordinates(self, coordinates: NDArray[np.float_]) -> None:
        self._coordinates = coordinates

    def __len__(self) -> int:
        return self._data.size

    def __getitem__(self, index: Union[int, slice]) -> NDArray[Any]:
        return self._data[index]

    def __iter__(self) -> Iterator[NDArray[Any]]:
        for i in range(len(self)):
            yield self[i]


class Residues:
    dtype = np.dtype({"names": ["resname", "resid"], "formats": ["<U10", "int"]})

    def __init__(self) -> None:
        self._data: NDArray[Any] = np.empty(0, dtype=self.dtype)

    @classmethod
    def from_gemmi(cls, gemmi_block: Dict[str, Any]) -> "Residues":
        cls_instance = cls.__new__(cls)

        # Create structured array
        resnames = np.array(gemmi_block["label_comp_id"])
        resids = np.array(gemmi_block["auth_seq_id"], dtype=int)

        data = np.empty(len(resnames), dtype=cls.dtype)
        data["resname"] = resnames
        data["resid"] = resids

        cls_instance._data = data

        return cls_instance

    @classmethod
    def from_mda(cls, mda_universe: UniverseType) -> "Residues":
        cls_instance = cls.__new__(cls)
        cls_instance._data = np.empty(len(mda_universe.atoms), dtype=cls_instance.dtype)
        cls_instance._data["resname"] = mda_universe.atoms.resnames
        cls_instance._data["resid"] = mda_universe.atoms.resids

        return cls_instance

    @property
    def resnames(self) -> NDArray[np.str_]:
        return self._data["resname"]

    @property
    def resids(self) -> NDArray[np.int32]:
        return self._data["resid"]

    def __len__(self) -> int:
        return len(self._data)

    def __getitem__(self, index: Union[int, slice]) -> NDArray[Any]:
        return self._data[index]

    def __iter__(self) -> Iterator[NDArray[Any]]:
        for i in range(len(self)):
            yield self[i]


class Chains:
    dtype = np.dtype(
        {"names": ["label", "auth", "id"], "formats": ["<U10", "<U10", "int"]}
    )

    def __init__(self, name: Optional[str] = None):
        self.name = name
        self._data: NDArray[Any] = np.empty(0, dtype=self.dtype)

        self.mapping: Dict[str, int] = {}

    @classmethod
    def from_gemmi(cls, gemmi_block: Dict[str, Any]) -> "Chains":
        cls_instance = cls.__new__(cls)

        # Create structured array
        labels: NDArray[np.str_] = np.array(gemmi_block["label_asym_id"])
        auths: NDArray[np.str_] = np.array(gemmi_block["auth_asym_id"])
        _, ids = np.unique(auths, return_inverse=True)  # type: ignore
        ids += 1

        data = np.empty(len(labels), dtype=cls.dtype)
        data["label"] = labels
        data["auth"] = auths
        data["id"] = ids

        cls_instance._data = data
        cls_instance.mapping = dict(zip(auths, ids))

        return cls_instance

    @classmethod
    def from_mda(cls, mda_universe: UniverseType) -> "Chains":
        cls_instance = cls.__new__(cls)
        cls_instance._data = np.empty(len(mda_universe.atoms), dtype=cls_instance.dtype)
        cls_instance._data["label"] = mda_universe.atoms.chainIDs
        cls_instance._data["auth"] = mda_universe.atoms.chainIDs
        _, cls_instance._data["id"] = np.unique(  # type: ignore
            cls_instance._data["auth"], return_inverse=True
        )
        cls_instance._data["id"] += 1

        cls_instance.mapping = dict(
            zip(cls_instance._data["auth"], cls_instance._data["id"])
        )

        return cls_instance

    @property
    def labels(self) -> NDArray[np.str_]:
        return self._data["label"]

    @property
    def auths(self) -> NDArray[np.str_]:
        return self._data["auth"]

    @property
    def ids(self) -> NDArray[np.int32]:
        return self._data["id"]

    def __len__(self) -> int:
        return len(self._data)

    def __getitem__(self, index: Union[int, slice]) -> NDArray[Any]:
        return self._data[index]

    def __iter__(self) -> Iterator[NDArray[Any]]:
        for i in range(len(self)):
            yield self[i]


class ARC:
    def __init__(
        self,
        obj: Union["GemmiLoader", "TopologyLoader"],
        site_data: Union[Dict[str, Any], AtomGroupType],
    ):
        obj_name: str = obj.__class__.__name__
        obj_map = self._obj_map(obj_name)
        self._atoms: Atoms = getattr(Atoms, obj_map)(site_data)
        self._residues: Residues = getattr(Residues, obj_map)(site_data)
        self._chains: Chains = getattr(Chains, obj_map)(site_data)

    def _obj_map(self, obj_name: str) -> str:
        mapping = {
            "GemmiLoader": "from_gemmi",
            "TopologyLoader": "from_mda",
        }
        return mapping[obj_name]

    def get_atom(self, index: int) -> "Atom":
        atom_info = self._atoms[index]
        residue_info = self._residues[index]
        chain_info = self._chains[index]

        atom_kwargs = {
            "name": atom_info["name"],
            "id": atom_info["id"],
            "element": atom_info["element"],
            "type": atom_info["type"],
            "resname": residue_info["resname"],
            "resid": residue_info["resid"],
            "chain_label": chain_info["label"],
            "chain_id": chain_info["id"],
        }

        return Atom(**atom_kwargs)

    @property
    def atoms(self) -> "Atoms":
        return self._atoms

    @property
    def residues(self) -> "Residues":
        return self._residues

    @property
    def chains(self) -> "Chains":
        return self._chains

    def __len__(self) -> int:
        return len(self._atoms)

    def __getitem__(self, index: Union[int, slice]) -> Union["Atom", List["Atom"]]:
        if isinstance(index, int):
            return self.get_atom(index)
        else:
            indices = np.arange(len(self))[index]  # type: ignore
            return [self.get_atom(i) for i in indices]

    def __iter__(self) -> Iterator["Atom"]:
        return (self.get_atom(i) for i in range(len(self)))


@dataclass
class Atom:
    # __slots__ = ["name", "id", "element", "type", "resname", "resid", "chain_label", "chain_id"]

    def __init__(self, **kwargs: Any):  # FIXME: specify the type of kwargs
        self.__dict__.update(kwargs)

    def __repr__(self) -> str:
        attrs = ", ".join(f"{k}={v}" for k, v in self.__dict__.items())
        return f"Atom({attrs})"
