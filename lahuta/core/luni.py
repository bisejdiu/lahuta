"""Defines the Luni class, which is the main class of the Lahuta package.
It is the entry point for all computations. It provides an interface for loading files,
or for initializing the Luni from an existing MDAnalysis.AtomGroup instance. This way
we can support both file-based loading and indicrectly support all MDAnalysis formats,
as well as provide support for reading MD trajectories.

Class:
    Luni: The main class of the Lahuta package.

Example:
-------
    universe = Luni(...)
    ns = universe.compute_neighbors()
    
"""

from typing import Any, Callable, List, Literal, Optional, Tuple, Union, overload

import MDAnalysis as mda
import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csc_array

from lahuta.config._atom_type_strings import BASE_AA_CONVERSION, RESIDUE_SYNONYMS
from lahuta.config.defaults import GEMMI_SUPPRTED_FORMATS
from lahuta.core._loaders import BaseLoader, GemmiLoader, TopologyLoader
from lahuta.core.arc import ARC
from lahuta.core.atom_assigner import AtomTypeAssigner
from lahuta.core.neighbor_finder import NeighborSearch
from lahuta.core.neighbors import NeighborPairs
from lahuta.core.topattrs import AtomAttrClassHandler  # This also imports VDWRadiiAtomAttr (which is needed)
from lahuta.lahuta_types.mdanalysis import AtomGroupType
from lahuta.lahuta_types.openbabel import MolType
from lahuta.utils.radii import v_radii_assignment

__all__ = ["LuniInputType", "Luni"]

LuniInputType = Union[AtomGroupType, str, List[str]]

class Luni:
    """The main class of the Lahuta package. It represents a universe of atoms and provides
    methods for computing various properties of the universe.

    The Luni class is the entry point for all computations. It provides an interface for loading files,
    or for initializing the Luni from an existing MDAnalysis.AtomGroup instance. This way we can support
    both file-based loading and indicrectly support all MDAnalysis formats, as well as provide support
    for reading MD trajectories.

    Args:
        *args (LuniInputType): Either an MDAnalysis.AtomGroup instance or a list of file names.

    Attributes:
        atom_types (csc_array): A sparse array containing the atom types of the universe.
        arc (ARC): The ARC instance used to load the files.
        sequence (str): The sequence of the universe.

    Raises:
        ValueError: If no input is provided or if invalid types of inputs are provided.
    """

    def __init__(self, *args: LuniInputType) -> None:
        self._mol: Optional[MolType] = None
        self._ready = False
        self._topattr_handler = AtomAttrClassHandler()
        self.atom_types: csc_array = csc_array((0, 0), dtype=np.int32)

        initializer = self._validate_input(*args)
        self._file_loader, self._mda = initializer(*args)

        assert self._mda is not None
        assert self._file_loader is not None

    def _validate_input(self, *args: LuniInputType) -> Callable[..., Tuple[BaseLoader, AtomGroupType]]:
        """Validate the input files for Luni initialization.

        Args:
        ----
        *args (LuniInputType): Either an MDAnalysis.AtomGroup instance or a list of file names.

        Raises:
            ValueError: If no input is provided or if invalid types of inputs are provided.

        Returns:
            func: A function to initialize the Luni either from an existing Luni or from files.
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
        """Initialize the universe from an existing Luni.

        Args:
            *args (LuniInputType): An MDAnalysis.AtomGroup instance.

        Returns:
            tuple: A tuple of the file loader and the AtomGroup instance.
        """
        _file_loader = TopologyLoader.from_mda(args[0])  # type: ignore
        _mda = _file_loader.to("mda")
        return _file_loader, _mda

    def _initialize_from_files(self, files: str) -> Tuple[BaseLoader, AtomGroupType]:
        """Initialize the universe from provided files.

        Args:
            files (str): The file name(s).

        Returns:
            tuple: A tuple of the file loader and the AtomGroup instance.
        """
        _file_loader = self._get_file_loader(files)
        _mda = _file_loader.to("mda")
        return _file_loader, _mda

    def _get_file_loader(self, *files: str) -> BaseLoader:
        """Retrieve the appropriate file loader.

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
            file_format, is_pdb = Luni.get_format(file_name)
            if file_format:
                return GemmiLoader(file_name, is_pdb=is_pdb)

        # If there are multiple files or the single file is not supported by GemmiLoader,
        # then use TopologyLoader
        return TopologyLoader(*files)

    def _extend_topology(self, attrname: str, values: NDArray[Any]) -> None:
        """Add new topology attributes to the Luni.

        Args:
            attrname (str): The name of the attribute.
            values (NDArray[Any]): The values of the attribute.
        """
        self._topattr_handler.init_topattr(attrname, attrname)
        self._mda.universe.add_TopologyAttr(attrname, values)

    def ready(self) -> None:
        """Prepare instance for computations by transforming the molecule and assigning atom types."""
        assert self._file_loader is not None
        self._mol = self._file_loader.to("mol")

        # TODO: remove array from the variable names by instead using type hints
        assert self.arc is not None
        atomtype_assigner = AtomTypeAssigner(self._mda, self._mol, legacy=False)
        ag_types = atomtype_assigner.assign_atom_types()
        og_atoms = self._mda.universe.atoms
        self.atom_types = ag_types

        # self._extend_topology("vdw_radii", v_radii_assignment(og_atoms.elements))
        self._mda.universe.add_TopologyAttr("vdw_radii", v_radii_assignment(og_atoms.elements))
        # for atom_type, value in AVAILABLE_ATOM_TYPES.items():
        #     self._extend_topology(atom_type.lower(), ag_types_array[:, value])

        self._ready = True

    # TODO: rename to find_neighbors
    def compute_neighbors(
        self,
        radius: float = 5.0,
        res_dif: int = 1,
    ) -> NeighborPairs:
        """Compute the neighbors of each atom in the Luni.

        This method calculates the neighbors for each atom based on the given radius and residue difference parameters.
        It returns an object of type `NeighborPairs` where each row in the underlying NumPy array contains the indices
        of the neighbors for the atom corresponding to that row.

        The method also ensures that the Luni instance is ready for computations by calling the `ready` method if needed.

        Args:
            radius (float, optional): The cutoff radius for considering two atoms as neighbors. Default is 5.0.
            res_dif (int, optional): The minimum difference in residue numbers for two atoms to be considered neighbors. Default is 1.

        Returns:
            NeighborPairs: An object containing a 2D NumPy array with shape (n_atoms, n_neighbors). Each row in the array
                        contains the indices of the neighbors for the atom corresponding to that row.
        """
        if not self._ready:
            self.ready()

        neighbors = NeighborSearch(self.to("mda"))
        pairs, distances = neighbors.compute(
            radius=radius,
            res_dif=res_dif,
        )

        return NeighborPairs(self.to("mda"), self.to("mol"), self.atom_types, pairs, distances)

    @property
    def sequence(self) -> str:
        """Retrieve the sequence of the Luni.

        This method retrieves the sequence of the Luni from the underlying MDA AtomGroup instance.
        It returns a NumPy array of shape (n_atoms,) containing the one-letter amino acid codes.

        Returns:
            NDArray[np.str_]: A NumPy array containing the one-letter amino acid codes of the Luni.
        """
        assert self.arc is not None
        three_letter_codes = self.to("mda").residues.resnames
        conversion_dict: dict[str, str] = {}
        for key, synonyms in RESIDUE_SYNONYMS.items():
            for synonym in synonyms:
                conversion_dict[synonym] = BASE_AA_CONVERSION[key]

        single_letter_codes = np.vectorize(lambda x: conversion_dict.get(x, 'X'))(three_letter_codes)

        return "".join(single_letter_codes)

    @staticmethod
    def get_format(file_name: str) -> Tuple[Union[str, None], bool]:
        """Retrieve the file format from a file name.

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
        """Convert the Luni to a different format.

        This method converts the internal representation of the Luni to the specified format.
        Currently, the supported formats are "mda" and "mol". If the desired format is "mol" and
        the Luni already has a "mol" representation, it returns the existing representation.
        Otherwise, it uses the `to` method of the file loader to perform the conversion.

        Args:
            fmt (str): The format to convert to. Currently supported formats are "mda" and "mol".

        Returns:
            Union[MolType, AtomGroupType]: A new Luni instance in the specified format.
        """
        if fmt not in {"mda", "mol"}:
            raise ValueError(f"Invalid format: {fmt}, must be one of 'mda' or 'mol'")

        if fmt == "mol" and self._mol is not None:
            return self._mol
        if fmt == "mol":
            self._mol = self._file_loader.to(fmt)
        return getattr(self, f"_{fmt}")  # type: ignore

    @property
    def arc(self) -> Union[None, ARC]:
        """Retrieve the ARC instance used to load the files.

        Returns:
            Union[None, ARC]: The ARC instance used to load the files.
        """
        return self._file_loader.arc

    def __repr__(self) -> str:
        return f"<Luni with {self._mda.n_atoms} atoms>"

    def __str__(self) -> str:
        return self.__repr__()
