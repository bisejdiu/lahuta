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

from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Optional, Type, Union, overload

import MDAnalysis as mda
import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csc_array, load_npz, save_npz
from typing_extensions import Self

from lahuta.config.atom_types import BASE_AA_CONVERSION, RESIDUE_SYNONYMS
from lahuta.config.atoms import PROT_ATOM_TYPES
from lahuta.config.defaults import GEMMI_SUPPRTED_FORMATS, MDA_SUPPORTED_FORMATS
from lahuta.core._loaders import BaseLoader, GemmiLoader, TopologyLoader
from lahuta.core.atom_assigner import AtomTypeAssigner
from lahuta.core.fn import GemmiNeighbors
from lahuta.core.neighbor_finder import NeighborSearch
from lahuta.core.neighbors import NeighborPairs
from lahuta.core.topattrs import AtomAttrClassHandler  # This also imports VDWRadiiAtomAttr (which is needed)
from lahuta.utils.array_utils import cross_interaction_indices

if TYPE_CHECKING:
    from lahuta.core.arc import ARC
    from lahuta.lahuta_types.mdanalysis import AtomGroupType, TrajectoryType
    from lahuta.lahuta_types.openbabel import MolType

__all__ = ["Luni"]

NeighborBackends: dict[str, Type[NeighborSearch] | Type[GemmiNeighbors]] = {
    "mda": NeighborSearch,
    "gemmi": GemmiNeighbors,
}


class Luni:
    """The main class of the Lahuta package. It represents a universe of atoms and provides
    methods for computing various properties of the universe.

    The Luni class is the entry point for all computations. It provides an interface for loading files,
    or for initializing the Luni from an existing MDAnalysis.AtomGroup instance. This way we can support
    both file-based loading and indicrectly support all MDAnalysis formats, as well as provide support
    for reading MD trajectories.

    Args:
        structure (str | AtomGroupType): The structure to load. Can be either a file path or an MDAnalysis AtomGroup.
        trajectories (Optional[str | list[str]]): Optional trajectories. Either a single path or a list of paths.

    Attributes:
        atom_types (csc_array): A sparse array containing the atom types of the universe.
        arc (ARC): The ARC instance used to load the files.
        sequence (str): The sequence of the universe.

    Raises:
        ValueError: If no input is provided or if invalid types of inputs are provided.
    """

    def __init__(
        self, structure: Union[str, Path, "AtomGroupType"], trajectories: Optional[str | list[str]] = None
    ) -> None:
        fmts: str | set[str] = ""
        self._file_loader: BaseLoader
        structure = str(structure) if isinstance(structure, Path) else structure
        match (structure, trajectories):
            case (mda.AtomGroup(atoms=s), None):
                self._file_loader = TopologyLoader.from_mda(s)  # type: ignore
            case (str(s), None):
                # Check if we can use GemmiLoader.
                file_format, is_pdb = Luni._check_gemmi_support(s)
                if file_format:
                    self._file_loader = GemmiLoader(s, is_pdb=is_pdb)
                elif s.upper().split(".")[-1] in MDA_SUPPORTED_FORMATS:
                    self._file_loader = TopologyLoader((s))
                else:
                    fmts = self._get_supported_fmts()
                    fmts = ", ".join(fmts)
                    raise ValueError(f"Unsupported format for structure: {s}! \nSupported formats are: {fmts}.")
            case (str(s), str(t)):
                # If trajectories are provided, use TopologyLoader
                self._file_loader = TopologyLoader(s, t)
            case _:
                fmts = self._get_supported_fmts()
                fmts = ", ".join(fmts)
                raise ValueError("Invalid input! \nSupported formats are: {fmts}.")

        self._mol: Optional["MolType"] = None
        self._mda = self._file_loader.to("mda")
        self.atom_types = csc_array((self._mda.universe.atoms.n_atoms, len(PROT_ATOM_TYPES)), dtype=np.int8)

        self._structure, self._trajectories = structure, trajectories

    @staticmethod
    def _get_supported_fmts() -> set[str]:
        return GEMMI_SUPPRTED_FORMATS.union({x.lower() for x in MDA_SUPPORTED_FORMATS})

    def _extend_topology(self, attrname: str, values: NDArray[Any]) -> None:
        """Add new topology attributes to the Luni.

        Args:
            attrname (str): The name of the attribute.
            values (NDArray[Any]): The values of the attribute.
        """
        topattr_handler = AtomAttrClassHandler()
        topattr_handler.init_topattr(attrname, attrname)
        self._mda.universe.add_TopologyAttr(attrname, values)

    def assing_atom_types(self) -> None:
        """Assign atom types to the Luni.

        This method assigns atom types to the Luni. It creates a sparse array of shape (n_atoms, n_atom_types)
        containing the atom types.

        Ensures that the Luni instance is ready for contact analysis.
        """
        if np.any(self.atom_types.data):
            return
        self._mol = self._file_loader.to("mol")
        atomtype_assigner = AtomTypeAssigner(self._mda, self._mol, legacy=False)
        self.atom_types = atomtype_assigner.assign_atom_types()

    def load_atom_types(self, filename: str, backend: Literal["scipy", "numpy"] = "numpy") -> None:
        """Load atom types from a file.

        This method loads the atom types from a file and stores them in the Luni object.

        Args:
            filename (str): The name of the file to load the atom types from.
            backend (Literal["scipy", "numpy"], optional): The backend to use for loading the atom types. \
            Default is "numpy".
        """
        self._mol = self._file_loader.to("mol")
        if backend == "scipy":
            self._load_using_scipy(filename)
        elif backend == "numpy":
            self._load_using_numpy(filename)
        else:
            raise ValueError(f"Invalid backend: {backend}, must be one of 'scipy' or 'numpy'")

    def unassign_atom_types(self) -> None:
        """Unassign atom types from the Luni.

        This method removes the atom types from the Luni.
        """
        self.atom_types *= 0

    def store_atom_types(
        self, filename: str = "sparse_matrix.npz", backend: Literal["scipy", "numpy"] = "scipy"
    ) -> None:
        """Store the atom types of the Luni to a file.

        This method stores the atom types of the Luni to a file.

        Args:
            filename (str): The name of the file to store the atom types in. Default is "sparse_matrix.npz".
            backend (Literal["scipy", "numpy"], optional): The backend to use for storing the atom types. \
            Default is "scipy".
        """
        assert np.any(self.atom_types.data), "Atom types have not been assigned yet!"
        if backend == "scipy":
            self._save_using_scipy(filename)
        elif backend == "numpy":
            self._save_using_numpy(filename)
        else:
            raise ValueError(f"Invalid backend: {backend}, must be one of 'scipy' or 'numpy'")

    def _save_using_scipy(self, filename: str) -> None:
        """Save the Luni to a file using SciPy.

        This method saves the Luni to a file using SciPy.

        Args:
            filename (str): The name of the file to save the Luni to.
        """
        save_npz(filename, self.atom_types)

    def _save_using_numpy(self, filename: str) -> None:
        """Save the Luni to a file using NumPy.

        This method saves the Luni to a file using NumPy.

        Args:
            filename (str): The name of the file to save the Luni to.
        """
        sparse_matrix = self.atom_types
        np.savez_compressed(
            filename,
            data=sparse_matrix.data,
            indices=sparse_matrix.indices,
            indptr=sparse_matrix.indptr,
            shape=sparse_matrix.shape,
        )

    def _load_using_scipy(self, filename: str) -> None:
        """Load the Luni from a file using SciPy.

        This method loads the Luni from a file using SciPy.

        Args:
            filename (str): The name of the file to load the Luni from.
        """
        self.atom_types = load_npz(filename)

    def _load_using_numpy(self, filename: str) -> None:
        """Load the Luni from a file using NumPy.

        This method loads the Luni from a file using NumPy.

        Args:
            filename (str): The name of the file to load the Luni from.
        """
        # loaded_matrix = csc_array((loaded['data'], loaded['indices'], loaded['indptr']), shape=loaded['shape'])
        sparse_matrix = np.load(filename)
        self.atom_types = csc_array((sparse_matrix["shape"]), dtype=np.int8)
        self.atom_types.data = sparse_matrix["data"]
        self.atom_types.indices = sparse_matrix["indices"]
        self.atom_types.indptr = sparse_matrix["indptr"]

    def compute_neighbors(  # noqa: PLR0913
        self,
        radius: float = 5.0,
        res_dif: int = 1,
        chain_type: Optional[Literal["inter", "intra"]] = None,
        image: Optional[Literal["inter", "intra"]] = None,
        target_spec: Optional[Self] = None,
        backend: Literal["mda", "gemmi"] = "mda",
        atom_types: bool = True,
    ) -> NeighborPairs:
        """Compute the neighbors of each atom in the Luni object.

        This method calculates the neighbors for each atom based on the given radius and residue difference parameters.
        It returns an object of type `NeighborPairs` where each row in the underlying NumPy array contains the indices
        of the neighbors for the atom corresponding to that row.

        Ensures that the Luni instance is ready for computations by calling the `ready` method if needed.

        Args:
            radius (float, optional): The cutoff radius for considering two atoms as neighbors. Default is 5.0.
            res_dif (int, optional): The minimum difference in residue numbers for two atoms to be \
            considered neighbors. Default is 1.
            chain_type (Literal["inter", "intra"], optional): The type of chain to keep. Default is None (keep all).
            image (Literal["inter", "intra"], optional): The type of image to keep. Default is None (keep all).
            backend (Literal["mda", "gemmi"], optional): The backend to use for computing neighbors. \
            Default is "mda".
            target_spec (Optional[Luni], optional): The target Luni object to compute neighbors with. \
            Default is None (compute neighbors within the same Luni object).
            atom_types (bool, optional): Whether to assign atom types before compute neighbors. Default is True.

        Returns:
            NeighborPairs: An object containing a 2D NumPy array with shape (n_atoms, n_neighbors). \
            Each row in the array contains the indices of the neighbors for the atom corresponding to that row.
        """
        if image is not None and backend != "gemmi":
            raise ValueError("Image filtering is only supported with the 'gemmi' backend!")

        self.assing_atom_types() if atom_types else self.unassign_atom_types()

        mda = self._mda
        if target_spec is not None:
            union_indices = np.union1d(self.indices, target_spec.indices)
            mda = self._mda.universe.atoms[union_indices]

        # neighbors = NeighborSearch(self.to("mda"))
        if backend == "gemmi":
            assert self._file_loader.structure is not None
            # TODO(bisejdiu): image is not being passed
            neighbors = GemmiNeighbors(mda, self._file_loader.structure)
        elif backend == "mda":
            neighbors = NeighborSearch(mda)
        else:
            raise ValueError(f"Invalid backend: {backend}, must be one of 'mda' or 'gemmi'")

        pairs, distances = neighbors.compute(
            radius=radius,
            res_dif=res_dif,
            chain_type=chain_type,
        )

        if target_spec is not None:
            cross_indices = cross_interaction_indices(pairs, self.indices, target_spec.indices)
            pairs, distances = pairs[cross_indices], distances[cross_indices]

        ns = NeighborPairs(self)
        ns.set_neighbors(pairs, distances)
        return ns

    # alias for compute_neighbors
    neighbors = compute_neighbors

    def filter(self, selection: str) -> Self:
        """Filter the Luni.

        This method filters the Luni based on the given selection string. It returns a new Luni instance
        containing the filtered atoms.

        Args:
            selection (str): The selection string to use for filtering.

        Returns:
            Luni: A new Luni instance containing the filtered atoms.
        """
        return self.__class__(self._mda.select_atoms(selection))

    def copy(self) -> Self:
        """Create a copy of this Luni instance.

        Returns:
            Luni: A copy of this Luni instance.
        """
        return self.__class__(self._mda)

    @property
    def sequence(self) -> str:
        """Retrieve the sequence of the Luni.

        This method retrieves the sequence of the Luni from the underlying MDA AtomGroup instance.
        It returns a NumPy array of shape (n_atoms,) containing the one-letter amino acid codes.

        Returns:
            NDArray[np.str_]: A NumPy array containing the one-letter amino acid codes of the Luni.
        """
        three_letter_codes = self.to("mda").select_atoms("protein").residues.resnames
        conversion_dict: dict[str, str] = {}
        for key, synonyms in RESIDUE_SYNONYMS.items():
            for synonym in synonyms:
                conversion_dict[synonym] = BASE_AA_CONVERSION[key]

        single_letter_codes = np.vectorize(lambda x: conversion_dict.get(x, "X"))(three_letter_codes)

        return "".join(single_letter_codes)

    @staticmethod
    def _check_gemmi_support(file_name: str) -> tuple[str | None, bool]:
        """Check if we can use Gemmi to load the file.

        Checks the file extension of the provided file name against the list of formats
        supported by GEMMI. If the extension matches a supported format,
        it returns the format and a boolean indicating whether the format is 'pdb' or 'cif'. If the file
        extension doesn't match any supported formats, it returns None and False.

        Args:
            file_name (str): The name of the file.

        Returns:
            tuple: A tuple containing the file format (str or None) and a boolean indicating if it is 'pdb' or 'cif'.
        """
        file_name_lower = file_name.lower()
        for fmt in GEMMI_SUPPRTED_FORMATS:
            if file_name_lower.endswith("." + fmt):
                is_pdb = fmt in {"pdb", "pdb.gz"}
                return fmt, is_pdb
        return None, False

    @overload
    def to(self, fmt: Literal["mda"]) -> "AtomGroupType":
        ...

    @overload
    def to(self, fmt: Literal["mol"]) -> "MolType":
        ...

    def to(self, fmt: Literal["mda", "mol"]) -> Union["MolType", "AtomGroupType"]:
        """Convert the Luni to a different format.

        This method converts the internal representation of the Luni to the specified format.
        Currently, the supported formats are "mda" and "mol". If the desired format is "mol" and
        the Luni already has a "mol" representation, it returns the existing representation.
        Otherwise, it uses the `to` method of the file loader to perform the conversion.

        Args:
            fmt (str): The format to convert to. Currently supported formats are "mda" and "mol".

        Returns:
            MolType | AtomGroupType: A new Luni instance in the specified format.
        """
        if fmt not in {"mda", "mol"}:
            raise ValueError(f"Invalid format: {fmt}, must be one of 'mda' or 'mol'")

        if fmt == "mol" and self._mol is not None:
            return self._mol
        if fmt == "mol":
            self._mol = self._file_loader.to(fmt)
        return getattr(self, f"_{fmt}")  # type: ignore

    @property
    def arc(self) -> "ARC":
        """Retrieve the ARC instance used to load the files.

        Returns:
            ARC: The ARC instance used to load the files.
        """
        assert self._file_loader.arc is not None, "Empty initialization is not currently supported!"
        return self._file_loader.arc

    @property
    def indices(self) -> NDArray[np.int32]:
        """Retrieve the indices of the atoms in the Luni object.

        Returns:
            NDArray[np.int32]: A NumPy array containing the indices of the atoms in the Luni object.
        """
        return self.arc.atoms.ids

    @property
    def ids(self) -> NDArray[np.int32]:
        """Retrieve the indices of the atoms in the Luni object.

        Returns:
            NDArray[np.int32]: A NumPy array containing the indices of the atoms in the Luni object.
        """
        return self.arc.atoms.ids

    @property
    def names(self) -> NDArray[np.str_]:
        """Retrieve the names of the atoms in the Luni object.

        Returns:
            NDArray[np.str_]: A NumPy array containing the names of the atoms in the Luni object.
        """
        return self.arc.atoms.names

    @property
    def elements(self) -> NDArray[np.str_]:
        """Retrieve the elements of the atoms in the Luni object.

        Returns:
            NDArray[np.str_]: A NumPy array containing the elements of the atoms in the Luni object.
        """
        return self.arc.atoms.elements

    @property
    def types(self) -> NDArray[np.str_]:
        """Retrieve the types of the atoms in the Luni object.

        Returns:
            NDArray[np.str_]: A NumPy array containing the types of the atoms in the Luni object.
        """
        return self.arc.atoms.types

    @property
    def coordinates(self) -> NDArray[np.float32]:
        """Retrieve the coordinates of the atoms in the Luni object.

        Returns:
            NDArray[np.float32]: A NumPy array containing the coordinates of the atoms in the Luni object.
        """
        return self.arc.atoms.coordinates

    @property
    def resnames(self) -> NDArray[np.str_]:
        """Retrieve the residue names of the atoms in the Luni object.

        Returns:
            NDArray[np.str_]: A NumPy array containing the residue names of the atoms in the Luni object.
        """
        return self.arc.residues.resnames

    @property
    def resids(self) -> NDArray[np.int32]:
        """Retrieve the residue IDs of the atoms in the Luni object.

        Returns:
            NDArray[np.int32]: A NumPy array containing the residue IDs of the atoms in the Luni object.
        """
        return self.arc.residues.resids

    @property
    def chainids(self) -> NDArray[np.int32]:
        """Retrieve the chain IDs of the atoms in the Luni object.

        Returns:
            NDArray[np.str_]: A NumPy array containing the chain IDs of the atoms in the Luni object.
        """
        return self.arc.chains.ids

    @property
    def chainlabels(self) -> NDArray[np.str_]:
        """Retrieve the chain labels of the atoms in the Luni object.

        Returns:
            NDArray[np.str_]: A NumPy array containing the chain labels of the atoms in the Luni object.
        """
        return self.arc.chains.labels

    @property
    def chainauths(self) -> NDArray[np.str_]:
        """Retrieve the chain auths of the atoms in the Luni object.

        Returns:
            NDArray[np.str_]: A NumPy array containing the chain auths of the atoms in the Luni object.
        """
        return self.arc.chains.auths

    @property
    def n_atoms(self) -> int:
        """Retrieve the number of atoms in the Luni object.

        Returns:
            int: The number of atoms in the Luni object.
        """
        return self.arc.atoms.n_atoms

    @property
    def n_residues(self) -> int:
        """Retrieve the number of residues in the Luni object.

        Returns:
            int: The number of residues in the Luni object.
        """
        return self._mda.residues.n_residues

    @property
    def n_chains(self) -> int:
        """Retrieve the number of chains in the Luni object.

        Returns:
            int: The number of chains in the Luni object.
        """
        return self._mda.n_segments

    @property
    def n_frames(self) -> int:
        """Retrieve the number of frames in the Luni object.

        Returns:
            int: The number of frames in the Luni object.
        """
        return self._mda.universe.trajectory.n_frames

    @property
    def trajectory(self) -> "TrajectoryType":
        """Retrieve the trajectory of the Luni object.

        Returns:
            Any: The trajectory of the Luni object.
        """
        return self._mda.universe.trajectory

    def __getstate__(self) -> dict[str, Any]:
        """Get state for pickling."""
        state = self.__dict__.copy()
        state["_mol"] = None
        return state

    def __setstate__(self, state: dict[str, Any]) -> None:
        """Set state for unpickling."""
        self.__dict__.update(state)

    def __reduce__(self) -> tuple[Type[Self], tuple[Any, ...], dict[str, Any]]:
        """Get state for pickling."""
        state = self.__dict__.copy()
        state["_mol"] = None
        return (self.__class__, (self._structure, self._trajectories), state)

    def __repr__(self) -> str:
        return f"<Luni with {self._mda.n_atoms} atoms>"

    def __str__(self) -> str:
        return self.__repr__()
