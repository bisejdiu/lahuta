"""Classes to model and manage the atomic-level structure of a biological system.
The atomic-level structure includes atoms, residues, and chains in the system, and their inter-relationships.

This module provides the following classes:

- `Atoms`: This class models and manages the properties of atoms in the system.
It provides functionalities to retrieve information about individual atoms,
as well as to iterate over and access data of multiple atoms.

- `Residues`: This class models and manages the properties of residues in the system.
Similar to the Atoms class, it provides functionalities to retrieve information about individual residues,
as well as to iterate over and access data of multiple residues.

- `Chains`: This class models and manages the properties of chains in the system.
Again, similar to the Atoms and Residues classes, it provides functionalities to retrieve information 
about individual chains, as well as to iterate over and access data of multiple chains.

- `ARC`: This class integrates the Atoms, Residues, and Chains (ARC) module, by creating and storing
instances of Atoms, Residues, and Chains classes. Depending on the type of the input object 
(either GemmiLoader or TopologyLoader), it uses the appropriate method from Atoms, Residues, and 
Chains classes to initialize them. It provides methods to retrieve atom, residue, 
and chain information individually or together, as well as methods to iterate over and access this data.

- `Atom`: This class models and manages the properties of an atom in the system. 
It is created with keyword arguments, which allows it to dynamically store any attributes passed.
The attributes include, but are not limited to, atom name, ID, element, type, 
residue name, residue ID, chain label, and chain ID.

Example:
    ``` py
    >>> from arc import ARC, GemmiLoader
    >>> loader = GemmiLoader("path_to_structure")
    >>> arc = ARC(loader, loader.atom_site_data)
    >>> atom_0 = arc.get_atom(0)
    >>> print(atom_0)
    Atom(name=..., id=..., element=..., type=..., resname=..., resid=..., chain_label=..., chain_id=...)
    ```

"""

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, Iterator, List, Optional, Union

import numpy as np
from numpy.typing import NDArray

from lahuta.lahuta_types.mdanalysis import AtomGroupType, UniverseType

if TYPE_CHECKING:
    from lahuta.core._loaders import GemmiLoader, TopologyLoader


class Atoms:
    """Atoms is a class that models and manages the properties and coordinates of atoms in a system.

    The Atoms class is part of the Atoms, Residues, Chains (ARC) module. It leverages structured
    numpy arrays to effectively handle, store and manipulate atomic information for speed optimization.
    The class stores atom-specific data like atom name, id, element type, and type, using the efficient
    and fast structured arrays provided by numpy.

    The class provides two class methods for initialization: `from_gemmi` and `from_mda`.
    These methods enable the creation of Atoms instances from Gemmi blocks or
    an MDAnalysis AtomGroup, respectively.

    The class also includes various properties to retrieve atomic attributes like names, types, ids,
    elements, and coordinates, and provides methods to iterate over and access the atomic data.

    Attributes:
        _data (NDArray[Any]): A numpy structured array storing atomic properties.
        _coordinates (NDArray[np.float32]): A numpy array storing the coordinates of atoms.

    Examples:
        ``` py
        >>> from arc import Atoms
        >>> atoms_from_gemmi = Atoms.from_gemmi(gemmi_block)
        >>> print(atoms_from_gemmi.names)
        ['C', 'O', 'N', ...]

        >>> atoms_from_mda = Atoms.from_mda(uv)
        >>> print(atoms_from_mda.ids)
        [1, 2, 3, ...]
        ```
    """

    dtype = np.dtype(
        {
            "names": ["name", "id", "element", "type"],
            "formats": ["<U10", "int", "<U10", "<U10"],
        }
    )

    def __init__(self) -> None:
        self._data: NDArray[Any] = np.empty(0, dtype=self.dtype)
        self._coordinates = np.zeros((0, 3), dtype=np.float32)

    @classmethod
    def from_gemmi(cls, gemmi_block: Dict[str, Any]) -> "Atoms":
        """Create an Atoms instance from a Gemmi block.

        Args:
            gemmi_block (Dict[str, Any]): A Gemmi block.

        Returns:
            Atoms: An Atoms instance.
        """
        cls_instance = cls.__new__(cls)

        # Create structured array
        label_atom_id: List[str] = gemmi_block["label_atom_id"]
        data = np.empty(len(label_atom_id), dtype=cls_instance.dtype)
        data["name"] = np.array(label_atom_id)
        data["id"] = np.arange(data["name"].size)
        data["element"] = np.array(gemmi_block.get("type_symbol"))
        data["type"] = np.array(gemmi_block.get("type_symbol"))

        cls_instance._data = data
        cls_instance._coordinates = np.zeros((0, 3), dtype=np.float32)

        return cls_instance

    @classmethod
    def from_mda(cls, uv: UniverseType) -> "Atoms":
        """Create an Atoms instance from an MDAnalysis Universe.

        Args:
            uv (UniverseType): An MDAnalysis Universe.

        Returns:
            Atoms: An Atoms instance.
        """
        cls_instance = cls.__new__(cls)
        data: NDArray[Any] = np.empty(len(uv.atoms), dtype=cls_instance.dtype)
        data["name"] = uv.atoms.names
        data["id"] = uv.atoms.ix
        data["element"] = uv.atoms.elements
        data["type"] = uv.atoms.types

        cls_instance._data = data
        cls_instance._coordinates = uv.atoms.positions

        return cls_instance

    @property
    def names(self) -> NDArray[np.str_]:
        """Atom names."""
        return self._data["name"]

    @property
    def types(self) -> NDArray[np.str_]:
        """Atom types."""
        return self._data["type"]

    @property
    def ids(self) -> NDArray[np.int32]:
        """Atom ids."""
        return self._data["id"]

    @property
    def elements(self) -> NDArray[np.str_]:
        """Atom elements."""
        return self._data["element"]

    @property
    def coordinates(self) -> NDArray[np.float32]:
        """3D coordinates of the atoms."""
        return self._coordinates

    @coordinates.setter
    def coordinates(self, coordinates: NDArray[np.float32]) -> None:
        self._coordinates = coordinates

    def __len__(self) -> int:
        return self._data.size

    def __getitem__(self, index: Union[int, slice]) -> NDArray[Any]:
        return self._data[index]

    def __iter__(self) -> Iterator[NDArray[Any]]:
        """Iterate over atoms."""
        for i in range(len(self)):
            yield self[i]


class Residues:
    """The Residues class models and manages the properties of residues in a system.

    The Residues class is part of the Atoms, Residues, Chains (ARC) module. It utilizes structured
    numpy arrays for effective handling, storage, and manipulation of residue-related data. This
    class specifically stores data related to residue names and residue IDs.

    Initialization of the Residues instances can be accomplished via two class methods:
    `from_gemmi` and `from_mda`. These methods enable the creation of Residues instances
    from Gemmi blocks or an MDAnalysis Universe respectively.

    The class includes properties to retrieve residue attributes like residue names and residue IDs,
    and provides methods to iterate over and access the residue data.

    Attributes:
        _data (NDArray[Any]): A numpy structured array storing residue properties.

    Examples:
        ``` py
        >>> from arc import Residues
        >>> residues_from_gemmi = Residues.from_gemmi(gemmi_block)
        >>> print(residues_from_gemmi.resnames)
        ['ALA', 'GLY', 'VAL', ...]

        >>> residues_from_mda = Residues.from_mda(uv)
        >>> print(residues_from_mda.resids)
        [1, 2, 3, ...]
        ```
    """

    dtype = np.dtype({"names": ["resname", "resid"], "formats": ["<U10", "int"]})

    def __init__(self) -> None:
        self._data: NDArray[Any] = np.empty(0, dtype=self.dtype)

    @classmethod
    def from_gemmi(cls, gemmi_block: Dict[str, Any]) -> "Residues":
        """Create a Residues instance from a Gemmi block.

        Args:
            gemmi_block (Dict[str, Any]): A Gemmi block.

        Returns:
            Residues: A Residues instance.
        """
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
        """Create a Residues instance from an MDAnalysis Universe.

        Args:
            mda_universe (UniverseType): An MDAnalysis Universe.

        Returns:
            Residues: A Residues instance.
        """
        cls_instance = cls.__new__(cls)
        cls_instance._data = np.empty(len(mda_universe.atoms), dtype=cls_instance.dtype)
        cls_instance._data["resname"] = mda_universe.atoms.resnames
        cls_instance._data["resid"] = mda_universe.atoms.resids

        return cls_instance

    @property
    def resnames(self) -> NDArray[np.str_]:
        """Residue names."""
        return self._data["resname"]

    @property
    def resids(self) -> NDArray[np.int32]:
        """Residue IDs."""
        return self._data["resid"]

    def __len__(self) -> int:
        return len(self._data)

    def __getitem__(self, index: Union[int, slice]) -> NDArray[Any]:
        return self._data[index]

    def __iter__(self) -> Iterator[NDArray[Any]]:
        """Iterate over residues."""
        for i in range(len(self)):
            yield self[i]


class Chains:
    """The Chains class models and manages the properties of chains in a system.

    The Chains class is part of the Atoms, Residues, Chains (ARC) module. It utilizes structured
    numpy arrays for effective handling, storage, and manipulation of chain-related data. This
    class specifically stores data related to chain labels, chain auth, and chain IDs.

    Initialization of the Chains instances can be accomplished via two class methods:
    `from_gemmi` and `from_mda`. These methods enable the creation of Chains instances from
    Gemmi blocks or an MDAnalysis Universe respectively.

    The class includes properties to retrieve chain attributes like labels, auths, and IDs,
    and provides methods to iterate over and access the chain data.

    Attributes:
        _data (NDArray[Any]): A numpy structured array storing chain properties.
        mapping (Dict[str, int]): Mapping from chain auth to chain IDs.

    Examples:
        ``` py
        >>> from arc import Chains
        >>> chains_from_gemmi = Chains.from_gemmi(gemmi_block)
        >>> print(chains_from_gemmi.labels)
        ['A', 'B', 'C', ...]

        >>> chains_from_mda = Chains.from_mda(uv)
        >>> print(chains_from_mda.auths)
        ['A', 'B', 'C', ...]
        ```
    """

    dtype = np.dtype({"names": ["label", "auth", "id"], "formats": ["<U10", "<U10", "int"]})

    def __init__(self, name: Optional[str] = None):
        self.name = name
        self._data: NDArray[Any] = np.empty(0, dtype=self.dtype)

        self.mapping: Dict[str, int] = {}

    @classmethod
    def from_gemmi(cls, gemmi_block: Dict[str, Any]) -> "Chains":
        """Create a Chains instance from a Gemmi block.

        Args:
            gemmi_block (Dict[str, Any]): A Gemmi block.

        Returns:
            Chains: A Chains instance.
        """
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
        """Create a Chains instance from an MDAnalysis Universe.

        Args:
            mda_universe (UniverseType): An MDAnalysis Universe.

        Returns:
            Chains: A Chains instance.
        """
        cls_instance = cls.__new__(cls)
        cls_instance._data = np.empty(len(mda_universe.atoms), dtype=cls_instance.dtype)
        cls_instance._data["label"] = mda_universe.atoms.chainIDs
        cls_instance._data["auth"] = mda_universe.atoms.chainIDs
        _, cls_instance._data["id"] = np.unique(cls_instance._data["auth"], return_inverse=True)  # type: ignore
        cls_instance._data["id"] += 1

        cls_instance.mapping = dict(zip(cls_instance._data["auth"], cls_instance._data["id"]))

        return cls_instance

    @property
    def labels(self) -> NDArray[np.str_]:
        """Chain labels."""
        return self._data["label"]

    @property
    def auths(self) -> NDArray[np.str_]:
        """Chain auths."""
        return self._data["auth"]

    @property
    def ids(self) -> NDArray[np.int32]:
        return self._data["id"]

    def __len__(self) -> int:
        return len(self._data)

    def __getitem__(self, index: Union[int, slice]) -> NDArray[Any]:
        return self._data[index]

    def __iter__(self) -> Iterator[NDArray[Any]]:
        """Iterate over chains."""
        for i in range(len(self)):
            yield self[i]


class ARC:
    """The ARC class models and manages the properties of a biological system consisting of atoms, residues, and chains.

    The ARC class integrates the Atoms, Residues, and Chains (ARC) module, by creating and storing
    instances of Atoms, Residues, and Chains classes.
    Depending on the type of the input object (either GemmiLoader or TopologyLoader), it uses the
    appropriate method from Atoms, Residues, and Chains classes to initialize them.

    It provides methods to retrieve atom, residue, and chain information individually or together,
    as well as methods to iterate over and access this data.

    Attributes:
        _atoms (Atoms): An instance of the Atoms class.
        _residues (Residues): An instance of the Residues class.
        _chains (Chains): An instance of the Chains class.

    Examples:
        ``` py
        >>> from arc import ARC, GemmiLoader
        >>> loader = GemmiLoader("path_to_structure")
        >>> arc = ARC(loader, loader.site_data)
        >>> atom_0 = arc.get_atom(0)
        >>> print(atom_0)
        Atom(name=..., id=..., element=..., type=..., resname=..., resid=..., chain_label=..., chain_id=...)

        >>> print(arc.atoms)
        Atoms(name=..., id=..., element=..., type=...)

        >>> print(arc.residues)
        Residues(resname=..., resid=...)
        ```
    """

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
        """Map the object name to the appropriate class method.

        Args:
            obj_name (str): The name of the object.

        Returns:
            str: The name of the class method.
        """
        mapping = {
            "GemmiLoader": "from_gemmi",
            "TopologyLoader": "from_mda",
        }
        return mapping[obj_name]

    def get_atom(self, index: int) -> "Atom":
        """Get an atom at a given index.

        Args:
            index (int): The index of the atom.

        Returns:
            Atom: An Atom instance.
        """
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
        """Atoms instance."""
        return self._atoms

    @property
    def residues(self) -> "Residues":
        """Residues instance."""
        return self._residues

    @property
    def chains(self) -> "Chains":
        """Chains instance."""
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
        """Iterate over atoms."""
        return (self.get_atom(i) for i in range(len(self)))


@dataclass
class Atom:
    """The Atom class models and manages the properties of an atom in a system.

    The Atom class stores information related to a single atom in a system.
    It's created using keyword arguments, which allows it to dynamically store any attributes passed.
    The attributes include, but are not limited to, atom name, ID, element, type,
    residue name, residue ID, chain label, and chain ID.

    The Atom class includes methods for representing and printing atom information.
    Note that all attributes are optional and can be passed as keyword arguments.

    Examples:
        ``` py
        >>> from arc import Atom
        >>> atom = Atom(name='CA', id=1, element='C', type='C.3', resname='ALA', resid=1, chain_label='A', chain_id=1)
        >>> print(atom)
        Atom(name=CA, id=1, element=C, type=C.3, resname=ALA, resid=1, chain_label=A, chain_id=1)
        ```
    """

    def __init__(self, **kwargs: Any):  # FIXME: specify the type of kwargs
        self.__dict__.update(kwargs)

    def __repr__(self) -> str:
        attrs = ", ".join(f"{k}={v}" for k, v in self.__dict__.items())
        return f"Atom({attrs})"
