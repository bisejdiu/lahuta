"""
Placeholder for the universe module.
"""

from typing import Any, Callable, List, Literal, Optional, Tuple, Union, overload

import MDAnalysis as mda
import numpy as np
from numpy.typing import NDArray

from lahuta.config.defaults import GEMMI_SUPPRTED_FORMATS
from lahuta.config.smarts import AVAILABLE_ATOM_TYPES
from lahuta.core._loaders import BaseLoader, GemmiLoader, TopologyLoader
from lahuta.core.arc import ARC
from lahuta.core.atom_assigner import AtomTypeAssigner
from lahuta.core.neighbor_finder import NeighborSearch
from lahuta.core.neighbors import NeighborPairs
from lahuta.core.topattrs import AtomAttrClassHandler
from lahuta.lahuta_types.mdanalysis import AtomGroupType
from lahuta.lahuta_types.openbabel import MolType
from lahuta.utils.radii import v_radii_assignment

LuniInputType = Union[AtomGroupType, str, List[str]]


class Universe:
    """
    This is the main class of the Lahuta package. It represents a universe of atoms and provides
    methods for computing various properties of the universe.

    Attributes:
        _mol (MolType, optional): A molecule in the universe.
        hbond_array (NoneType): Array to store hydrogen bonds.
        _ready (bool): State of the universe, whether it's ready for computations.
        _mapping (NDArray[np.int64]): Maps atom indices to their positions in a flat, 1D array.
        _topattr_handler (AtomAttrClassHandler): Handles atom attributes.
        _file_loader (BaseLoader, optional): Handles file loading.
        _mdag (AtomGroupType, optional): Represents a group of atoms in the Universe.

    Methods:
        _validate_input(*args: LuniInputType): Validates the input files.
        _initialize_from_universe(*args: LuniInputType): Initializes the universe from existing Universe.
        _initialize_from_files(files: str): Initializes the universe from provided files.
        _get_file_loader(*files: str): Retrieves the appropriate file loader.
        _extend_topology(attrname: str, values: NDArray[Any]): Adds new topology attributes to the Universe.
        _build_atom_mapping(ag: AtomGroupType): Builds a mapping of atom indices.
        ready(): Prepares instance for computations by transforming the molecule and assigning atom types.
        compute_neighbors(radius: float, res_dif: int): Computes the neighbors of each atom in the Universe.
        get_format(file_name: str): Retrieves the file format from a file name.
        to(fmt: Literal["mda", "mol"]): Converts the Universe to a different format.
        arc(): Retrieves the ARC object from the file loader.

    """

    def __init__(self, *args: LuniInputType) -> None:
        self._mol: Optional[MolType] = None
        self.hbond_array = None
        self._ready = False
        self._mapping: NDArray[np.int64] = np.array([], dtype=np.int64)
        self._topattr_handler = AtomAttrClassHandler()
        # self._file_loader: Optional[BaseLoader] = None
        # self._mdag: Optional[AtomGroupType] = None

        initializer = self._validate_input(*args)
        self._file_loader, self._mdag = initializer(*args)

        assert self._mdag is not None
        assert self._file_loader is not None

    def _validate_input(self, *args: LuniInputType) -> Callable[..., Tuple[BaseLoader, AtomGroupType]]:
        """
        Validates the input files for Universe initialization.

        Args:
            *args (LuniInputType): Either an MDAnalysis.AtomGroup instance or a list of file names.

        Raises:
            ValueError: If no input is provided or if invalid types of inputs are provided.

        Returns:
            func: A function to initialize the Universe either from an existing Universe or from files.
        """

        if not args:
            raise ValueError("No input provided")

        if isinstance(args[0], mda.AtomGroup):
            if len(args) != 1:
                raise ValueError("When passing an MDAnalysis.AtomGroup instance, no other arguments are allowed")
            return self._initialize_from_universe

        if not all(isinstance(arg, str) for arg in args):
            raise ValueError("Input must be either an MDAnalysis.AtomGroup instance or a list of file names")
        return self._initialize_from_files

    def _initialize_from_universe(self, *args: LuniInputType) -> Tuple[BaseLoader, AtomGroupType]:
        """
        Initializes the universe from an existing Universe.

        Args:
            *args (LuniInputType): An MDAnalysis.AtomGroup instance.

        Returns:
            tuple: A tuple of the file loader and the AtomGroup instance.
        """

        _file_loader = TopologyLoader.from_mda(args[0])  # type: ignore
        _mdag = _file_loader.to("mda")
        return _file_loader, _mdag

    def _initialize_from_files(self, files: str) -> Tuple[BaseLoader, AtomGroupType]:
        """
        Initializes the universe from provided files.

        Args:
            files (str): The file name(s).

        Returns:
            tuple: A tuple of the file loader and the AtomGroup instance.
        """

        _file_loader = self._get_file_loader(files)
        _mdag = _file_loader.to("mda")
        return _file_loader, _mdag

    def _get_file_loader(self, *files: str) -> BaseLoader:
        """
        Retrieves the appropriate file loader.

        Args:
            *files (str): The file name(s).

        Raises:
            ValueError: If no file name is provided or if multiple files are provided.

        Returns:
            BaseLoader: The appropriate file loader.
        """

        # GemmiLoader can only handle one file and its format should be supported
        if len(files) == 1:
            file_name = files[0]
            file_format, is_pdb = Universe.get_format(file_name)
            if file_format:
                return GemmiLoader(file_name, is_pdb=is_pdb)

        # If there are multiple files or the single file is not supported by GemmiLoader,
        # then use TopologyLoader
        return TopologyLoader(*files)

    def _extend_topology(self, attrname: str, values: NDArray[Any]) -> None:
        """
        Adds new topology attributes to the Universe.

        Args:
            attrname (str): The name of the attribute.
            values (NDArray[Any]): The values of the attribute.
        """

        self._topattr_handler.init_topattr(attrname, attrname)
        self._mdag.universe.add_TopologyAttr(attrname, values)

    # def select_atoms(self, *args, **kwargs) -> mda.AtomGroup:
    #     return self.atoms.select_atoms(*args, **kwargs)

    def _build_atom_mapping(self, ag: AtomGroupType) -> NDArray[np.int64]:
        """
        Builds a mapping of atom indices.

        Args:
            ag (AtomGroupType): The AtomGroup instance.

        Returns:
            NDArray[np.int64]: The atom mapping.
        """

        max_index = np.max(ag.universe.atoms.indices)
        atom_mapping = np.full(max_index + 1, -1, dtype=np.int64)
        atom_mapping[ag.indices] = np.arange(ag.n_atoms)
        return atom_mapping

    def ready(self) -> None:
        """
        Prepares instance for computations by transforming the molecule and assigning atom types.

        Raises:
            ValueError: If the Universe is already ready.
        """

        assert self._file_loader is not None
        self._mol = self._file_loader.to("mol")
        self._mapping = self._build_atom_mapping(self.to("mda").universe.atoms)

        # TODO: remove array from the variable names by instead using type hints
        atomtype_assigner = AtomTypeAssigner(self._mdag, self._mol, self._mapping)
        ag_types = atomtype_assigner.assign_atom_types()
        og_atoms = self._mdag.universe.atoms

        reference_array = np.zeros((og_atoms.n_atoms, ag_types.shape[1]))

        ix = self._mdag.indices
        full_ag_atypes = reference_array.copy()
        full_ag_atypes[ix] = ag_types

        self._extend_topology("vdw_radii", v_radii_assignment(og_atoms.elements))
        for atom_type, value in AVAILABLE_ATOM_TYPES.items():
            self._extend_topology(atom_type.lower(), full_ag_atypes[:, value])

        self._ready = True

    # TODO: rename to find_neighbors
    def compute_neighbors(
        self,
        radius: float = 5.0,
        res_dif: int = 1,
    ) -> NeighborPairs:
        """
        Compute the neighbors of each atom in the Universe.

        This method calculates the neighbors for each atom based on the given radius and residue difference parameters.
        It returns an object of type `NeighborPairs` where each row in the underlying NumPy array contains the indices
        of the neighbors for the atom corresponding to that row.

        The method also ensures that the Universe instance is ready for computations by calling the `ready` method if needed.

        Args:
            radius (float, optional): The cutoff radius for considering two atoms as neighbors. Default is 5.0.
            res_dif (int, optional): The minimum difference in residue numbers for two atoms to be considered neighbors. Default is 1.

        Returns:
            NeighborPairs: An object containing a 2D NumPy array with shape (n_atoms, n_neighbors). Each row in the array
                        contains the indices of the neighbors for the atom corresponding to that row.

        Raises:
            AssertionError: If the Universe instance is not ready for computations.

        Note:
            The actual computation of the neighbors is delegated to an instance of the `NeighborSearch` class.
        """

        if not self._ready:
            self.ready()

        neighbors = NeighborSearch(self.to("mda"))
        pairs, distances = neighbors.compute(
            radius=radius,
            res_dif=res_dif,
        )

        return NeighborPairs(self.to("mda"), self.to("mol"), pairs, distances)

    @staticmethod
    def get_format(file_name: str) -> Tuple[Union[str, None], bool]:
        """
        Retrieve the file format from a file name.

        This static method checks the file extension of the provided file name against the list of formats
        supported by GEMMI (stored in `GEMMI_SUPPRTED_FORMATS`). If the extension matches a supported format,
        it returns the format and a boolean indicating whether the format is 'pdb' or 'pdb.gz'. If the file
        extension doesn't match any supported formats, it returns None and False.

        Args:
            file_name (str): The name of the file.

        Returns:
            tuple: A tuple containing the file format (str or None) and a boolean indicating if it is 'pdb' or 'pdb.gz'.
        """
        file_name_lower = file_name.lower()
        for fmt in GEMMI_SUPPRTED_FORMATS:
            if file_name_lower.endswith("." + fmt):
                is_pdb = fmt in {"pdb", "pdb.gz"}
                return fmt, is_pdb
        return None, False

    @overload
    def to(self, fmt: Literal["mda"]) -> AtomGroupType:
        ...

    @overload
    def to(self, fmt: Literal["mol"]) -> MolType:
        ...

    def to(self, fmt: Literal["mda", "mol"]) -> Union[MolType, AtomGroupType]:
        """
        Convert the Universe to a different format.

        This method converts the internal representation of the Universe to the specified format.
        Currently, the supported formats are "mda" and "mol". If the desired format is "mol" and
        the Universe already has a "mol" representation, it returns the existing representation.
        Otherwise, it uses the `to` method of the file loader to perform the conversion.

        Args:
            fmt (str): The format to convert to. Currently supported formats are "mda" and "mol".

        Returns:
            Union[MolType, AtomGroupType]: A new Universe instance in the specified format.
        """
        if fmt == "mol" and self._mol is not None:
            return self._mol
        return self._file_loader.to(fmt)

    @property
    def arc(self) -> Union[None, ARC]:
        return self._file_loader.arc

    def __repr__(self) -> str:
        return f"<Lahuta Universe with {self._mdag.n_atoms} atoms>"

    def __str__(self) -> str:
        return self.__repr__()
