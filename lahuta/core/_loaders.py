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


from abc import ABC, abstractmethod
from typing import Any, Iterator, Literal, Optional, overload

import gemmi
import MDAnalysis as mda
import numpy as np
import pandas as pd
from numpy.typing import NDArray

from lahuta.core.arc import ARC, Atoms, Chains, Residues
from lahuta.core.obmol import OBMol
from lahuta.lahuta_types.mdanalysis import AtomGroupType, UniverseType
from lahuta.lahuta_types.openbabel import MolType
from lahuta.utils.radii import v_radii_assignment


class BaseLoader(ABC):
    """Abstract base class providing a blueprint for loading and handling biological structure data.

    `BaseLoader` forms the basis for other loader classes in the lahuta library, dedicated
    to handling different formats of biological structural data. It includes the initialization
    of atomic, residue, and chain data as well as their transformation into different formats.

    The class takes a file path as an input and reads the biological structure data. It offers
    various properties to access the data, as well as methods to convert the loaded data into
    different output formats. It also provides an iterator to iterate over atoms, residues, and chains.

    Args:
        file_path (str): The file path from which to load the biological structure data.

    Attributes:
        file_path (str): The file path from which the biological structure data is loaded.
        _chains (Optional[Chains]): An instance of Chains class containing chain data.
        _residues (Optional[Residues]): An instance of Residues class containing residue data.
        _atoms (Optional[Atoms]): An instance of Atoms class containing atom data.
        structure: A structure object containing biological structure data, currently initialized to None.
        arc (Optional[ARC]): An instance of ARC class providing combined access to atoms, residues, and chains data.
    """

    def __init__(self, file_path: str):
        self.file_path = file_path
        self._chains = None
        self._residues = None
        self._atoms = None

        self.structure = None
        self.arc: Optional[ARC] = None

    @property
    def n_atoms(self) -> int:
        """The number of atoms in the loaded biological structure data."""
        return len(self.atoms)

    @property
    def chains(self) -> Chains:
        """The chains in the loaded biological structure data."""
        if self.arc is None:
            raise ValueError("arc has not been initialized")
        return self.arc.chains

    @property
    def residues(self) -> Residues:
        """The residues in the loaded biological structure data."""
        if self.arc is None:
            raise ValueError("arc has not been initialized")
        return self.arc.residues

    @property
    def atoms(self) -> Atoms:
        """The atoms in the loaded biological structure data."""
        if self.arc is None:
            raise ValueError("arc has not been initialized")
        return self.arc.atoms

    @overload
    def to(self, fmt: Literal["mda"]) -> AtomGroupType:
        ...

    @overload
    def to(self, fmt: Literal["mol"]) -> MolType:
        ...

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

    def __iter__(self) -> Iterator[Atoms | Residues | Chains]:
        """Iterate over atoms, residues, and chains."""
        yield self.atoms
        yield self.residues
        yield self.chains


class GemmiLoader(BaseLoader):
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
        super().__init__(file_path)
        if is_pdb:
            structure: Any = gemmi.read_pdb(self.file_path)  # type: ignore
            block: Any = structure.make_mmcif_document().sole_block()
        else:
            block: Any = gemmi.cif.read(self.file_path).sole_block()  # type: ignore
            structure: Any = gemmi.make_structure_from_block(block)  # type: ignore

        self.structure = structure
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
        coords_array: NDArray[np.float32] = np.zeros((self.n_atoms, 3), dtype=np.float32)
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

        resnames, resids, chain_ids = (uniques["resnames"], uniques["resids"], uniques["chain_ids"])

        # Create a new Universe
        uv: UniverseType = mda.Universe.empty(
            n_atoms=self.arc.atoms.ids.size,
            n_residues=uniques.size,
            n_segments=chain_ids.size,
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
        uv.add_TopologyAttr("segids", chain_ids)
        uv.add_TopologyAttr("ids", self.arc.atoms.ids)

        uv.atoms.positions = self.arc.atoms.coordinates

        return uv.atoms

    def to_mda(self) -> AtomGroupType:
        """Convert the loaded biological structure data into an MDAnalysis AtomGroup object."""
        return self.ag

    def to_mol(self) -> MolType:
        """Convert the loaded biological structure data into an OpenBabel Mol object."""
        assert self.arc is not None, "arc has not been initialized"
        assert self.structure is not None, "structure should not be None"
        obmol = OBMol()
        obmol.create_mol(
            self.arc,
            self.structure.connections,
        )
        assert obmol.mol is not None
        return obmol.mol


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
                trajectories = (trajectories)
            self.ag.universe.load_new(trajectories, format=None, in_memory=False)

        self.ag.universe.add_TopologyAttr("vdw_radii", v_radii_assignment(universe.atoms.elements))
        self.arc = ARC(self, self.ag)  # positions are set when using mda.Universe

    def to_mda(self) -> AtomGroupType:
        """Convert the loaded biological structure data into an MDAnalysis AtomGroup object."""
        return self.ag

    def to_mol(self) -> MolType:
        """Convert the loaded biological structure data into an OpenBabel Mol object."""
        assert self.arc is not None, "arc has not been initialized"
        obmol = OBMol()
        obmol.create_mol(
            self.arc,
            self.structure,
        )
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
        top_loader.ag = ag.copy()
        top_loader.ag._u = ag.universe.copy()  # noqa: SLF001
        top_loader.structure = None

        top_loader.arc = ARC(top_loader, top_loader.ag)

        return top_loader
