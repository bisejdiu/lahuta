from typing import TYPE_CHECKING, Union

import numpy as np
from MDAnalysis.core.universe import Universe

if TYPE_CHECKING:
    from lahuta.core.loaders import GemmiLoader, TopologyLoader


class Residues:
    def __init__(self, name=None):
        self.name = name

        self._resnames = np.array([], dtype="<U10")
        self._resids = np.array([], dtype=int)
        self._resindices = np.array([], dtype=int)
        self._icodes = np.array([], dtype="<U10")

    @classmethod
    def from_gemmi(cls, gemmi_block):
        cls_instance = cls.__new__(cls)
        cls_instance._resnames = np.array(gemmi_block.get("label_comp_id"))
        cls_instance._resids = np.array(gemmi_block.get("auth_seq_id"), dtype=int)
        cls_instance._icodes = np.array(
            gemmi_block.get("pdbx_PDB_ins_code"), dtype=bool
        )
        _, cls_instance._resindices = np.unique(
            cls_instance._resids, return_inverse=True
        )

        return cls_instance

    @classmethod
    def from_mda(cls, mda_universe):
        cls_instance = cls.__new__(cls)
        cls_instance._resnames = mda_universe.atoms.resnames
        cls_instance._resids = mda_universe.atoms.resids
        cls_instance._resindices = mda_universe.atoms.resindices
        cls_instance._icodes = mda_universe.atoms.icodes

        return cls_instance

    @property
    def resnames(self):
        return self._resnames

    @property
    def resids(self):
        return self._resids

    @property
    def resindices(self):
        return self._resindices

    @property
    def icodes(self):
        return self._icodes

    def __len__(self):
        return len(self._resnames)

    def __getitem__(self, index):
        return self._resnames[index], self._resids[index], self._resindices[index]

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


class Atoms:
    dtype = np.dtype(
        {
            "names": ["name", "id", "element", "type"],
            "formats": ["<U10", "int", "<U10", "<U10"],
        }
    )

    def __init__(self, name=None):
        self.name = name

        # Define dtype for an atom
        # self.dtype = np.dtype(
        #     {
        #         "names": ["name", "id", "element", "type"],
        #         "formats": ["<U10", "int", "<U10", "<U10"],
        #     }
        # )

        self._data = np.empty(0, dtype=self.dtype)

    @classmethod
    def from_gemmi(cls, gemmi_block):
        cls_instance = cls.__new__(cls)

        # Create structured array
        label_atom_id = gemmi_block.get("label_atom_id")
        data = np.empty(len(label_atom_id), dtype=cls_instance.dtype)
        data["name"] = np.array(label_atom_id)
        data["id"] = np.array(gemmi_block.get("id"), dtype=int)
        data["element"] = np.array(gemmi_block.get("type_symbol"))
        data["type"] = np.array(gemmi_block.get("type_symbol"))

        cls_instance._data = data

        return cls_instance

    @classmethod
    def from_mda(cls, mda_universe):
        cls_instance = cls.__new__(cls)
        cls_instance._data = np.empty(len(mda_universe.atoms), dtype=cls_instance.dtype)
        cls_instance._data["name"] = mda_universe.atoms.names
        cls_instance._data["id"] = mda_universe.atoms.ids
        cls_instance._data["element"] = mda_universe.atoms.elements
        cls_instance._data["type"] = mda_universe.atoms.types

        return cls_instance

    @property
    def names(self):
        return self._data["name"]

    @property
    def types(self):
        return self._data["type"]

    @property
    def ids(self):
        return self._data["id"]

    @property
    def elements(self):
        return self._data["element"]

    def __len__(self):
        return len(self._data)

    def __getitem__(self, index):
        return self._data[index]

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


# class Atoms:
#     def __init__(self, name=None):
#         self.name = name

#         self._names = np.array([], dtype="<U10")
#         self._ids = np.array([], dtype=int)
#         self._ix = np.array([], dtype=int)  # re-indexed ids
#         self._elements = np.array([], dtype="<U10")
#         self._types = np.array([], dtype="<U10")

#     @classmethod
#     def from_gemmi(cls, gemmi_block):
#         cls_instance = cls.__new__(cls)
#         cls_instance._names = np.array(gemmi_block.get("label_atom_id"))
#         cls_instance._ids = np.array(gemmi_block.get("id"), dtype=int)
#         cls_instance._ix = cls_instance._ids.copy()
#         cls_instance._elements = np.array(gemmi_block.get("type_symbol"))
#         cls_instance._types = np.array(gemmi_block.get("type_symbol"))

#         return cls_instance

#     @classmethod
#     def from_mda(cls, mda_universe):
#         cls_instance = cls.__new__(cls)
#         cls_instance._names = mda_universe.atoms.names
#         cls_instance._ids = mda_universe.atoms.ids
#         cls_instance._ix = np.arange(cls_instance._ids.size)
#         cls_instance._elements = mda_universe.atoms.elements
#         cls_instance._types = mda_universe.atoms.types

#         return cls_instance

#     @property
#     def names(self):
#         return self._names

#     @property
#     def types(self):
#         return self._types

#     @property
#     def ids(self):
#         return self._ids

#     @property
#     def ix(self):
#         return self._ix

#     @property
#     def elements(self):
#         return self._elements

#     def __len__(self):
#         return len(self._names)

#     def __getitem__(self, index):
#         return self._names[index], self._ids[index], self._elements[index]

#     def __iter__(self):
#         for i in range(len(self)):
#             yield self[i]


class Chains:
    def __init__(self, name=None):
        self.name = name

        self._labels = np.array([], dtype="<U10")
        self._auths = np.array([], dtype="<U10")
        self._ids = np.array([], dtype=int)

        self.mapping = {}

    @classmethod
    def from_gemmi(cls, gemmi_block):
        cls_instance = cls.__new__(cls)
        cls_instance._labels = np.array(gemmi_block.get("label_asym_id"))
        cls_instance._auths = np.array(gemmi_block.get("auth_asym_id"))
        _, cls_instance._ids = np.unique(cls_instance._auths, return_inverse=True)
        cls_instance._ids += 1

        cls_instance.mapping = dict(zip(cls_instance._auths, cls_instance._ids))

        return cls_instance

    @classmethod
    def from_mda(cls, mda_universe):
        cls_instance = cls.__new__(cls)
        cls_instance._labels = mda_universe.atoms.chainIDs
        cls_instance._auths = mda_universe.atoms.chainIDs
        _, cls_instance._ids = np.unique(cls_instance._auths, return_inverse=True)
        cls_instance._ids += 1

        cls_instance.mapping = dict(zip(cls_instance._auths, cls_instance._ids))

        return cls_instance

    @property
    def labels(self):
        return self._labels

    @property
    def auths(self):
        return self._auths

    @property
    def ids(self):
        return self._ids

    def __len__(self):
        return len(self._labels)

    def __getitem__(self, index):
        return self._labels[index], self._ids[index]

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


class CRA:
    def __init__(self, obj: Union["GemmiLoader", "TopologyLoader"], site_data):
        obj_name = obj.__class__.__name__
        obj_map = self._obj_map(obj_name)
        print("obj_map", obj_map)
        self._atoms = getattr(Atoms, obj_map)(site_data)
        self._residues = getattr(Residues, obj_map)(site_data)
        self._chains = getattr(Chains, obj_map)(site_data)

    def _obj_map(self, obj_name):
        mapping = {"GemmiLoader": "from_gemmi", "TopologyLoader": "from_mda"}
        return mapping[obj_name]

    @property
    def atoms(self):
        return self._atoms

    @property
    def residues(self):
        return self._residues

    @property
    def chains(self):
        return self._chains

    def __iter__(self):
        return zip(self._chains, self._residues, self._atoms)
