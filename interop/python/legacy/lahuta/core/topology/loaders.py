"""Classes to load and manage biological structure data.

The two main classes are `GemmiLoader` and `TopologyLoader`, each designed to work with a
specific library, Gemmi and MDAnalysis respectively, to load and handle structure data.

Both loaders extend from an abstract base class `BaseLoader` which outlines the necessary
methods and properties all loader classes should have.

The module also provides a unified way to manage atoms, residues, and chains information
through the ARC class instance. It includes various methods for querying the data and converting
them into other common bioinformatics and chemoinformatics formats.

Classes:
    ```
    BaseLoader: Abstract base class that provides a template for the loader classes.
    GemmiLoader: Class to load and manage biological structure data using the Gemmi library.
    TopologyLoader: Class to load and manage biological structure data using the MDAnalysis library.
    ```
"""

from abc import ABC, abstractmethod, abstractproperty
from typing import Any, Iterator, Literal, NamedTuple, Optional, Protocol, overload

import gemmi
import MDAnalysis as mda
import numpy as np

from lahuta.lib._lahuta import IR, LahutaSystem
import pandas as pd
from numpy.typing import NDArray

from lahuta.lib import cLuni
from lahuta._types.mdanalysis import AtomGroupType, UniverseType
from lahuta._types.openbabel import MolType
from lahuta.utils.radii import v_radii_assignment

from .arc import ARC, Atom, Atoms, Chains, Residues
from .obmol import OBMol


class PartnerData(NamedTuple):
    """Store gemmi.ConnectionList atom partner data."""

    atom_name: str
    chain_name: str
    res_id_seq_num: int
    res_id_name: str


class ConnectionData(NamedTuple):
    """Store gemmi.ConnectionList atom connection data."""

    partner1: PartnerData
    partner2: PartnerData


# FIX: load_file should take the loader as an argument
def load_file(file_path: str) -> LahutaSystem | AtomGroupType:
    # return LahutaSystemLoader(file_path) if file_path.endswith(".luni") else TopologyLoader(file_path)
    lahuta_supported = {".pdb", ".cif", ".pdb.gz", ".cif.gz"}
    if file_path.endswith(tuple(lahuta_supported)):
        return LahutaSystemLoader(file_path).luni
    return TopologyLoader(file_path).ag


# class Loader(Protocol):
#     def __init__(self):
#         self.input: LahutaSystem
#
#     @staticmethod
#     def to_ir(input: k


class Loader(Protocol):
    def __init__(self):
        self.luni: LahutaSystem

    # @staticmethod
    # def from_ir(ir: IR) -> "Loader": ...
    #
    # @staticmethod
    # def to_ir() -> IR: ...

    @property
    @abstractmethod
    def n_atoms(self) -> int:
        """The number of atoms in the loaded biological structure data."""
        ...

    @property
    @abstractmethod
    def indices(self) -> NDArray[np.int32]:
        """The number of atoms in the loaded biological structure data."""
        ...

    @property
    @abstractmethod
    def names(self) -> NDArray[np.str_]:
        """The number of atoms in the loaded biological structure data."""
        ...

    @property
    @abstractmethod
    def elements(self) -> NDArray[np.str_]:
        """The number of atoms in the loaded biological structure data."""
        ...

    @property
    @abstractmethod
    def positions(self) -> NDArray[np.float32]:
        """The number of atoms in the loaded biological structure data."""
        ...

    @property
    @abstractmethod
    def resnames(self) -> NDArray[np.str_]:
        """The residue names in the loaded biological structure data."""
        ...

    @property
    @abstractmethod
    def resids(self) -> NDArray[np.int32]:
        """The residue IDs in the loaded biological structure data."""
        ...


# Intermediate Representation
# class IRType:
#     def __init__(self, data: Any):
#         self.data = data


class BaseLoader(ABC):
    """Abstract base class providing a blueprint for loading and handling biological structure data.

    `BaseLoader` forms the basis for other loader classes in the lahuta library, dedicated
    to handling different formats of biological structural data. It includes the initialization
    of atomic, residue, and chain data as well as their transformation into different formats.

    The class takes a file path as an input and reads the biological structure data. It offers
    various properties to access the data, as well as methods to convert the loaded data into
    different output formats. It also provides an iterator to iterate over atoms, residues, and chains.
    """

    def __init__(self, file_path: str):
        self.file_path = file_path
        self._data: LahutaSystem | AtomGroupType
        self.luni: LahutaSystem

    @property
    @abstractmethod
    def n_atoms(self) -> int:
        """The number of atoms in the loaded biological structure data."""
        ...

    @property
    @abstractmethod
    def indices(self) -> NDArray[np.int32]:
        """The number of atoms in the loaded biological structure data."""
        ...

    @property
    @abstractmethod
    def ids(self) -> NDArray[np.int32]:
        """The number of atoms in the loaded biological structure data."""
        ...

    @property
    @abstractmethod
    def names(self) -> NDArray[np.str_]:
        """The number of atoms in the loaded biological structure data."""
        ...

    @property
    @abstractmethod
    def elements(self) -> NDArray[np.str_]:
        """The number of atoms in the loaded biological structure data."""
        ...

    @property
    @abstractmethod
    def coordinates(self) -> NDArray[np.float32]:
        """The number of atoms in the loaded biological structure data."""
        ...

    # @property
    # @abstractmethod
    # def chains(self) -> Chains:
    #     """The chains in the loaded biological structure data."""
    #     ...

    @property
    @abstractmethod
    def resnames(self) -> NDArray[np.str_]:
        """The residue names in the loaded biological structure data."""
        ...

    @property
    @abstractmethod
    def resids(self) -> NDArray[np.int32]:
        """The residue IDs in the loaded biological structure data."""
        ...

    @overload
    def to(self, fmt: Literal["mda"]) -> AtomGroupType: ...

    @overload
    def to(self, fmt: Literal["mol"]) -> MolType: ...

    def to(self, fmt: Literal["mol", "mda"]) -> MolType | AtomGroupType:
        """Convert the loaded biological structure data into a different format.

        Args:
            fmt (str): The format to which the biological structure data should be converted.

        Returns:
            MolType | AtomGroupType: The biological structure data in the specified format.

        Raises:
            ValueError: If the specified format is not supported.
        """
        method_str = f"to_{fmt}"
        if hasattr(self, method_str):
            return getattr(self, method_str)()  # type: ignore

        raise ValueError(f"Object type {fmt} is not supported")

    @abstractmethod
    def to_mda(self) -> AtomGroupType:
        """Convert the loaded biological structure data into an MDAnalysis AtomGroup object."""

    @abstractmethod
    def to_mol(self) -> MolType:
        """Convert the loaded biological structure data into an OpenBabel Mol object."""


class LahutaSystemLoader(BaseLoader):
    """Class for loading biological structure data using the Lahuta C++ library."""

    def __init__(self, file_path: str):
        self.file_path = file_path
        self.luni: LahutaSystem = cLuni(file_path)
        print("Done with C++")

    @property
    def n_atoms(self) -> int:
        return self.luni.n_atoms

    @property
    def ids(self) -> NDArray[np.int32]:
        return self.luni.indices

    @property
    def indices(self) -> NDArray[np.int32]:
        return self.luni.indices

    @property
    def names(self) -> NDArray[np.str_]:
        return self.luni.names

    @property
    def elements(self) -> NDArray[np.str_]:
        return self.luni.elements

    @property
    def resnames(self) -> NDArray[np.str_]:
        return self.luni.resnames

    @property
    def resids(self) -> NDArray[np.int32]:
        return self.luni.resids

    @property
    def resindices(self) -> NDArray[np.int32]:
        return self.luni.resindices

    @property
    def chainlabels(self) -> NDArray[np.int32]:
        return self.luni.chainlabels

    @property
    def coordinates(self) -> NDArray[np.float32]:
        return self.luni.coordinates()

    def _create_mda(self) -> AtomGroupType:
        """Create an MDAnalysis AtomGroup object from the loaded data.

        Returns:
            AtomGroupType: An MDAnalysis AtomGroup object.
        """
        # Create a structured array to ensure unique values for each combination of resname, resid, and chain_id
        # import time

        # start = time.time()
        struct_arr = np.rec.fromarrays(  # type: ignore
            [self.resnames, self.resids, self.chainlabels],
            names=str("resnames, resids, chain_ids"),  # type: ignore
        )
        # print("Time taken to create structured array", time.time() - start)
        # start = time.time()

        # Use factorize to get the labels and unique values
        resindices, uniques = pd.factorize(struct_arr)
        # print("Time taken to factorize", time.time() - start)
        # start = time.time()
        resnames, resids, chain_ids = uniques["resnames"], uniques["resids"], uniques["chain_ids"]
        # print("Time taken to get uniques", time.time() - start)
        # start = time.time()
        # print("resindices", resindices)
        # print("resindices", resindices.shape, self.n_atoms)
        # print("chain_labels", self.chainlabels, self.chainlabels.shape, np.unique(self.chainlabels))
        # print("chain_labels", chain_ids, chain_ids.shape)

        _, int_array = np.unique(chain_ids, return_inverse=True)
        # print("Time taken to get unique and inverse", time.time() - start)
        # start = time.time()

        # Convert 0-based labels to 1-based labels
        int_array += 1
        # print("chain_ids", int_array, int_array.shape)

        # Create a new Universe
        uv: UniverseType = mda.Universe.empty(
            n_atoms=self.n_atoms,
            n_residues=uniques.size,
            n_segments=self.chainlabels.size,
            atom_resindex=resindices,
            residue_segindex=int_array,
            trajectory=True,
        )
        # print("Time taken to create Universe", time.time() - start)
        # start = time.time()

        # Add topology attributes
        uv.add_TopologyAttr("names", self.names)
        uv.add_TopologyAttr("type", self.names)
        uv.add_TopologyAttr("elements", self.elements)
        uv.add_TopologyAttr("vdw_radii", v_radii_assignment(self.elements))
        # uv.add_TopologyAttr("vdw_radii", v_radii)
        uv.add_TopologyAttr("resnames", resnames)
        uv.add_TopologyAttr("resids", resids)
        uv.add_TopologyAttr("chainIDs", self.chainlabels)
        uv.add_TopologyAttr("ids", self.indices)
        # uv.add_TopologyAttr("tempfactors", self.arc.atoms.b_isos)
        # print("Time taken to add topology attributes", time.time() - start)
        # start = time.time()

        uv.atoms.positions = self.coordinates
        uv.filename = self.file_path
        # print("Time taken to set positions", time.time() - start)

        return uv.atoms

    def to_mda(self) -> AtomGroupType:
        """Convert the loaded biological structure data into an MDAnalysis AtomGroup object."""
        return self._create_mda()

    def to_mol(self) -> Any:
        return None


class GemmiLoader:
    """Class for loading biological structure data using the gemmi library.

    `GemmiLoader` is a subclass of `BaseLoader` and provides the implementation to read,
    convert, and manage biological structure data specifically from the gemmi library.
    The class takes a file path as an input and reads the biological structure data either
    in PDB format or in the format specified in the file. It also offers methods to convert
    the loaded data into different output formats.

    Args:
        file_path (str): The file path from which to load the biological structure data.
        is_pdb (bool): A flag to indicate whether the file is in PDB format. Defaults to False.

    Attributes:
        structure: A gemmi structure object containing biological structure data.
        arc (Optional[ARC]): An instance of ARC class providing combined access to atoms,
                             residues, and chains data.
        ag (AtomGroupType): An MDAnalysis AtomGroup object created from the loaded data.
    """

    def __init__(self, file_path: str, is_pdb: bool = False):
        # super().__init__(file_path)
        self.file_path = file_path
        if is_pdb:
            structure = gemmi.read_pdb(self.file_path)
            block = structure.make_mmcif_document().sole_block()
        else:
            block = gemmi.cif.read(self.file_path).sole_block()
            structure = gemmi.make_structure_from_block(block)

        self.conn_store = self.store_connections(structure.connections)
        atom_site_data: dict[str, Any] = block.get_mmcif_category("_atom_site.")

        self.arc = ARC(self, atom_site_data)
        self.arc.atoms.coordinates = self.extract_positions(atom_site_data)

        self.ag: AtomGroupType = self._create_mda()

    def extract_positions(self, atom_site_data: dict[str, Any]) -> NDArray[np.float32]:
        """Extract the coordinates from the atom_site data.

        Args:
            atom_site_data (dict[str, Any]): The atom_site data from the gemmi structure object.

        Returns:
            NDArray[np.float32]: The coordinates of the atoms in the structure.

        Raises:
            ValueError: If the atom_site data does not contain the required information.
        """
        # coords_array: NDArray[np.float32] = np.zeros((self.n_atoms, 3), dtype=np.float32)
        coords_array: NDArray[np.float32] = np.zeros((len(atom_site_data.get("id")), 3), dtype=np.float32)
        coords_array[:, 0] = atom_site_data.get("Cartn_x")
        coords_array[:, 1] = atom_site_data.get("Cartn_y")
        coords_array[:, 2] = atom_site_data.get("Cartn_z")

        return coords_array

    def _create_mda(self) -> AtomGroupType:
        """Create an MDAnalysis AtomGroup object from the loaded data.

        Returns:
            AtomGroupType: An MDAnalysis AtomGroup object.
        """
        # Create a structured array to ensure unique values for each combination of resname, resid, and chain_id
        assert self.arc is not None, "arc has not been initialized"
        struct_arr = np.rec.fromarrays(  # type: ignore
            [self.arc.residues.resnames, self.arc.residues.resids, self.arc.chains.ids],
            names=str("resnames, resids, chain_ids"),
        )

        # Use factorize to get the labels and unique values
        resindices, uniques = pd.factorize(struct_arr)
        resnames, resids, chain_ids = uniques["resnames"], uniques["resids"], uniques["chain_ids"]

        # Create a new Universe
        uv: UniverseType = mda.Universe.empty(
            n_atoms=self.arc.atoms.ids.size,
            n_residues=uniques.size,
            n_segments=self.arc.chains.ids.size,
            atom_resindex=resindices,
            residue_segindex=chain_ids,
            trajectory=True,
        )

        # Add topology attributes
        uv.add_TopologyAttr("names", self.arc.atoms.names)
        uv.add_TopologyAttr("type", self.arc.atoms.types)
        uv.add_TopologyAttr("elements", self.arc.atoms.elements)
        uv.add_TopologyAttr("vdw_radii", v_radii_assignment(self.arc.atoms.elements))
        uv.add_TopologyAttr("resnames", resnames)
        uv.add_TopologyAttr("resids", resids)
        uv.add_TopologyAttr("chainIDs", self.arc.chains.auths)
        uv.add_TopologyAttr("ids", self.arc.atoms.ids)
        uv.add_TopologyAttr("tempfactors", self.arc.atoms.b_isos)

        uv.atoms.positions = self.arc.atoms.coordinates
        uv.filename = self.file_path

        return uv.atoms

    def to_mda(self) -> AtomGroupType:
        """Convert the loaded biological structure data into an MDAnalysis AtomGroup object."""
        return self.ag

    def to_mol(self) -> MolType:
        """Convert the loaded biological structure data into an OpenBabel Mol object."""
        assert self.arc is not None, "arc has not been initialized"
        obmol = OBMol()
        obmol.create_mol(
            self.arc,
            self.conn_store,
        )
        assert obmol.mol is not None
        return obmol.mol

    def store_connections(self, connections: gemmi.ConnectionList) -> list[ConnectionData]:
        """Store the atom connections from the gemmi structure object.

        Args:
            connections (gemmi.ConnectionList): The atom connections from the gemmi structure object.

        Returns:
            list[ConnectionData]: A list of ConnectionData objects containing the atom connections.
        """
        conn_store = []
        for conn in connections:
            prt1 = PartnerData(
                atom_name=conn.partner1.atom_name,
                chain_name=conn.partner1.chain_name,
                res_id_seq_num=conn.partner1.res_id.seqid.num,  # type: ignore
                res_id_name=conn.partner1.res_id.name,
            )
            prt2 = PartnerData(
                atom_name=conn.partner2.atom_name,
                chain_name=conn.partner2.chain_name,
                res_id_seq_num=conn.partner2.res_id.seqid.num,  # type: ignore
                res_id_name=conn.partner2.res_id.name,
            )
            conn_store.append(ConnectionData(prt1, prt2))

        return conn_store


class TopologyLoader(BaseLoader):
    """Class for loading and managing biological structure data using the MDAnalysis library.

    `TopologyLoader` is a subclass of `BaseLoader` and provides the implementation to read,
    convert, and manage biological structure data specifically from the MDAnalysis library.
    This class accepts one or more file paths and loads the structure data into an
    MDAnalysis Universe. It also offers methods to convert the loaded data into different output formats.

    Args:
        *paths (str): One or more file paths from which to load the biological structure data.

    Attributes:
        ag (AtomGroupType): An MDAnalysis AtomGroup object created from the loaded data.
        arc (Optional[ARC]): An instance of ARC class providing combined access to atoms,
                             residues, and chains data.

    """

    def __init__(self, structure: str, *trajectories: str):
        file_path: str = structure
        super().__init__(file_path)
        universe = mda.Universe(self.file_path)
        self.ag: AtomGroupType = universe.atoms
        assert self.ag is not None
        if trajectories:
            if isinstance(trajectories, str):
                trajectories = trajectories
            self.ag.universe.load_new(trajectories, format=None, in_memory=False)

        self.ag.universe.add_TopologyAttr("vdw_radii", v_radii_assignment(universe.atoms.elements))
        self.ag.universe.add_TopologyAttr("element", np.char.upper(universe.atoms.elements.astype(str)))
        self.ag.universe.filename = self.file_path
        # self.ag.elements = np.char.capitalize(self.ag.elements)
        self.arc = ARC(self, self.ag)  # positions are set when using mda.Universe

    def to_mda(self) -> AtomGroupType:
        """Convert the loaded biological structure data into an MDAnalysis AtomGroup object."""
        return self.ag

    def to_mol(self) -> MolType:
        """Convert the loaded biological structure data into an OpenBabel Mol object."""
        assert self.arc is not None, "arc has not been initialized"
        obmol = OBMol()
        obmol.create_mol(self.arc)
        assert obmol.mol is not None
        return obmol.mol

    @classmethod
    def from_mda(cls, ag: AtomGroupType) -> "TopologyLoader":
        """Create a new instance of the TopologyLoader class from an MDAnalysis AtomGroup object.

        Args:
            ag (AtomGroupType): An MDAnalysis AtomGroup object.

        Returns:
            TopologyLoader: A new instance of the TopologyLoader class with the AtomGroup data copied.
        """
        top_loader = cls.__new__(cls)
        ag.universe.add_TopologyAttr("vdw_radii", v_radii_assignment(ag.universe.atoms.elements))
        ag.universe.add_TopologyAttr("element", np.char.upper(ag.universe.atoms.elements.astype(str)))

        top_loader.ag = ag.copy()
        top_loader.ag._u = ag.universe.copy()  # noqa: SLF001
        top_loader.ag.universe.filename = ag.universe.filename

        top_loader.arc = ARC(top_loader, top_loader.ag)

        return top_loader

    @property
    def n_atoms(self) -> int:
        return self.ag.n_atoms

    @property
    def indices(self) -> NDArray[np.int32]:
        return self.ag.indices

    @property
    def ids(self) -> NDArray[np.int32]:
        return self.ag.ids

    @property
    def resnames(self) -> NDArray[np.str_]:
        return self.ag.resnames

    @property
    def resids(self) -> NDArray[np.int32]:
        return self.ag.resids

    @property
    def names(self) -> NDArray[np.str_]:
        return self.ag.names

    @property
    def elements(self) -> NDArray[np.str_]:
        return self.ag.elements

    @property
    def coordinates(self) -> NDArray[np.float32]:
        return self.ag.positions
