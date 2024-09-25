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

import time
from lahuta.core.topology.loaders import GemmiLoader
from lahuta.lib import cLuni, cAtomType

from lahuta.config.atom_types import BASE_AA_CONVERSION, RESIDUE_SYNONYMS
from lahuta.config.atoms import PROT_ATOM_TYPES
from lahuta.config.defaults import GEMMI_SUPPRTED_FORMATS, MDA_SUPPORTED_FORMATS
from lahuta.core.assigners import AtomTypeAssigner
from lahuta.core.neighbors import BaseNeighborSearch, NeighborPairs
from lahuta.core.neighbors.backends import GemmiNeighborSearch, MDAnalysisNeighborSearch
from lahuta.core.topology import LahutaCPPLoader, TopologyLoader
from lahuta.utils.array_utils import cross_interaction_indices

if TYPE_CHECKING:
    from lahuta._types.mdanalysis import AtomGroupType, TrajectoryType
    from lahuta._types.openbabel import MolType
    from lahuta.core.topology import BaseLoader
    from lahuta.core.topology.arc import ARC

__all__ = ["Luni"]

NeighborBackends: dict[str, Type[BaseNeighborSearch]] = {
    "mda": MDAnalysisNeighborSearch,
    "gemmi": GemmiNeighborSearch,
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
        self,
        structure: "str | Path | AtomGroupType",
        trajectories: str | list[str] | None = None,
        b_iso_name: str = "tempfactor",
    ) -> None:
        fmts: str | set[str] = ""
        self._file_loader: BaseLoader
        self._input_structure: str | Path | "AtomGroupType" = structure
        structure = str(structure) if isinstance(structure, Path) else structure
        match (structure, trajectories):
            case (mda.AtomGroup(atoms=s), None):
                self._file_loader = TopologyLoader.from_mda(s)  # type: ignore
            case (str(s), None):
                # Check if we can use GemmiLoader.
                file_format, is_pdb = Luni._check_gemmi_support(s)
                if file_format:
                    self._file_loader = LahutaCPPLoader(s)
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
                raise ValueError(
                    f"Invalid input! {structure=} and {trajectories=} are not valid inputs. \nSupported formats are: {fmts}."
                )

        self._mol: Optional["MolType"] = None
        self._mda = self._file_loader.to("mda")
        self.atom_types = csc_array((self._mda.universe.atoms.n_atoms, len(PROT_ATOM_TYPES)), dtype=np.int8)

        # self._luni = cLuni(self._file_loader.file_path)
        # self._at = self._luni.get_atom_types()

        self._structure, self._trajectories = structure, trajectories

        self.b_iso_name = b_iso_name
        if self.b_iso_name != "tempfactor":
            self.extend_topology(self.b_iso_name, self._mda.universe.atoms.tempfactors)

    @staticmethod
    def _get_supported_fmts() -> set[str]:
        return GEMMI_SUPPRTED_FORMATS.union({x.lower() for x in MDA_SUPPORTED_FORMATS})

    def extend_topology(self, attrname: str, values: NDArray[Any]) -> None:
        """Add new topology attributes to the Luni.

        Args:
            attrname (str): The name of the attribute.
            values (NDArray[Any]): The values of the attribute.
        """
        from .topology.topattrs import AtomAttrClassHandler

        topattr_handler = AtomAttrClassHandler()
        topattr_handler.init_topattr(attrname, attrname)
        self._mda.universe.add_TopologyAttr(attrname, values)

    def assing_atom_types(self) -> None:
        # pass

        #     """Assign atom types to the Luni.
        #
        #     This method assigns atom types to the Luni. It creates a sparse array of shape (n_atoms, n_atom_types)
        #     containing the atom types.
        #
        #     Ensures that the Luni instance is ready for contact analysis.
        #     """
        #     self.atom_types = csc_array((self._mda.universe.atoms.n_atoms, len(PROT_ATOM_TYPES)), dtype=np.int8)
        #     atom_types = self._file_loader.luni.get_atom_types()
        #     print("Atom types: ", atom_types)
        #     for i, at in enumerate(atom_types):
        #         self.atom_types[i, at] = 1

        # if np.any(self.atom_types.data):
        #     return
        tmp_loader = GemmiLoader(self._file_loader.file_path)
        self._mol = tmp_loader.to_mol()
        atomtype_assigner = AtomTypeAssigner(self._mda, self._mol, legacy=False)
        self.atom_types = atomtype_assigner.assign_atom_types()
        self._mda = tmp_loader.to_mda()

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

        self.assing_atom_types()  # if atom_types else self.unassign_atom_types()

        mda = self._mda
        if target_spec is not None:
            union_indices = np.union1d(self.indices, target_spec.indices)
            mda = self._mda.universe.atoms[union_indices]

        # neighbors = MDAnalysisNeighborSearch(mda)

        neighbors = self._file_loader.luni._find_neighbors(radius)
        neighbors = neighbors.remove_adjascent_pairs(res_dif)
        # print("ff", ff.get_pairs().shape)
        pairs, distances = neighbors.get_pairs(), np.array(neighbors.get_distances_sq())
        # pairs, distances = neighbors.get_pairs(), neighbors.get_distances()

        # print("------>", p.shape, d.shape)

        # pairs, distances = neighbors.compute(
        #     radius=radius,
        #     res_dif=res_dif,
        #     chain_type=chain_type,
        # )
        # print("compare with: ", pairs.shape, distances.shape)

        if target_spec is not None:
            cross_indices = cross_interaction_indices(pairs, self.indices, target_spec.indices)
            pairs, distances = pairs[cross_indices], distances[cross_indices]

        ns = NeighborPairs(self)
        ns.set_neighbors(pairs, distances)
        ns._ns = neighbors
        return ns

    # alias for compute_neighbors
    neighbors = compute_neighbors

    def filter(self, selection: str) -> Self:
        """Filter atoms from this Luni using a selection string to a new Luni.

        `lahuta` uses the MDAnalysis selection language to filter atoms. It defines
        several additional keywords in addition to the standard MDAnalysis keywords.

        **Dynamic selections**
            Not tested & not documented & probably not working.

        Args:
            selection (str): The selection string to use for filtering.

        Returns:
            Luni: A new Luni instance containing the filtered atoms.
        """
        if self.b_iso_name != "tempfactor":
            selection = selection.replace(self.b_iso_name, "tempfactor")
        return self.__class__(self._mda.select_atoms(selection))

    def remove_water(self) -> Self:
        """Remove water molecules from the Luni.

        This method removes water molecules from the Luni. It returns a new Luni instance
        containing the filtered atoms.

        Returns:
            Luni: A new Luni instance containing the filtered atoms.
        """
        return self.filter("not water")

    def remove_ions(self) -> Self:
        """Remove ions from the Luni.

        This method removes ions from the Luni. It returns a new Luni instance
        containing the filtered atoms.

        Returns:
            Luni: A new Luni instance containing the filtered atoms.
        """
        return self.filter("not ion")

    def remove_ligands(self) -> Self:
        """Remove ligands from the Luni.

        This method removes ligands from the Luni. It returns a new Luni instance
        containing the filtered atoms.

        Returns:
            Luni: A new Luni instance containing the filtered atoms.
        """
        return self.filter("not ligand")

    def remove_sugars(self) -> Self:
        """Remove sugars from the Luni.

        This method removes sugars from the Luni. It returns a new Luni instance
        containing the filtered atoms.

        Returns:
            Luni: A new Luni instance containing the filtered atoms.
        """
        return self.filter("not sugar")

    def remove_lipids(self) -> Self:
        """Remove lipids from the Luni.

        This method removes lipids from the Luni. It returns a new Luni instance
        containing the filtered atoms.

        Returns:
            Luni: A new Luni instance containing the filtered atoms.
        """
        return self.filter("not lipid")

    def remove_nucleic(self) -> Self:
        """Remove nucleic acids from the Luni.

        This method removes nucleic acids from the Luni. It returns a new Luni instance
        containing the filtered atoms.

        Returns:
            Luni: A new Luni instance containing the filtered atoms.
        """
        return self.filter("not nucleic")

    def remove_protein(self) -> Self:
        """Remove protein from the Luni.

        This method removes protein from the Luni. It returns a new Luni instance
        containing the filtered atoms.

        Returns:
            Luni: A new Luni instance containing the filtered atoms.
        """
        return self.filter("not protein")

    def copy(self) -> Self:
        """Create a copy of this Luni instance.

        Returns:
            Luni: A copy of this Luni instance.
        """
        return self.__class__(self._mda)

    def deepcopy(self) -> Self:
        """Create a deep copy of this Luni instance.

        Returns:
            Luni: A deep copy of this Luni instance.
        """
        mda_copy = self._mda.copy()
        _, resindices = np.unique(mda_copy.resindices, return_inverse=True)

        uv = mda.Universe.empty(
            n_atoms=mda_copy.n_atoms,
            n_residues=mda_copy.n_residues,
            n_segments=mda_copy.n_segments,
            atom_resindex=resindices,
            residue_segindex=mda_copy.residues.segindices,
            trajectory=True,
        )

        uv.add_TopologyAttr("names", mda_copy.names)
        uv.add_TopologyAttr("types", mda_copy.types)
        uv.add_TopologyAttr("elements", mda_copy.elements)
        uv.add_TopologyAttr("resnames", mda_copy.residues.resnames)
        uv.add_TopologyAttr("resids", mda_copy.residues.resids)
        uv.add_TopologyAttr("chainIDs", mda_copy.chainIDs)
        uv.add_TopologyAttr("ids", mda_copy.ids)
        uv.add_TopologyAttr("tempfactors", mda_copy.tempfactors)

        assert uv.atoms is not None
        uv.atoms.positions = mda_copy.positions
        uv.filename = mda_copy.universe.filename

        return self.__class__(uv.atoms)

    def sequence(self, use_synonyms: bool = False) -> str:
        """Retrieve the sequence of the Luni.

        This method retrieves the sequence of the Luni from the underlying MDA AtomGroup instance.
        It returns a NumPy array of shape (n_atoms,) containing the one-letter amino acid codes.

        Returns:
            NDArray[np.str_]: A NumPy array containing the one-letter amino acid codes of the Luni.
        """
        res_dict: dict[str, str] = {}
        three_letter_codes = self.to("mda").select_atoms("protein").residues.resnames
        if use_synonyms:
            for key, synonyms in RESIDUE_SYNONYMS.items():
                for synonym in synonyms:
                    res_dict[synonym] = BASE_AA_CONVERSION[key]

        res_dict = res_dict if use_synonyms else BASE_AA_CONVERSION
        unknown_res_name = "X" if use_synonyms else ""
        single_letter_codes = np.vectorize(lambda x: res_dict.get(x, unknown_res_name))(three_letter_codes)

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
    def to(self, fmt: Literal["mda"]) -> "AtomGroupType": ...

    @overload
    def to(self, fmt: Literal["mol"]) -> "MolType": ...

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

    # @property
    # def arc(self) -> "ARC":
    #     """Retrieve the ARC instance used to load the files.
    #
    #     Returns:
    #         ARC: The ARC instance used to load the files.
    #     """
    #     assert self._file_loader.arc is not None, "Empty initialization is not currently supported!"
    #     return self._file_loader.arc

    # @property
    # def indices(self) -> NDArray[np.int32]:
    #     """Retrieve the indices of the atoms in the Luni object.
    #
    #     Returns:
    #         NDArray[np.int32]: A NumPy array containing the indices of the atoms in the Luni object.
    #     """
    #     return self.arc.atoms.ids

    @property
    def indices(self) -> NDArray[np.int32]:
        """Retrieve the indices of the atoms in the Luni object.

        Returns:
            NDArray[np.int32]: A NumPy array containing the indices of the atoms in the Luni object.
        """
        return self._file_loader.indices

    # @property
    # def ids(self) -> NDArray[np.int32]:
    #     """Retrieve the indices of the atoms in the Luni object.
    #
    #     Returns:
    #         NDArray[np.int32]: A NumPy array containing the indices of the atoms in the Luni object.
    #     """
    #     return self.arc.atoms.ids

    @property
    def ids(self) -> NDArray[np.int32]:
        """Retrieve the indices of the atoms in the Luni object.

        Returns:
            NDArray[np.int32]: A NumPy array containing the indices of the atoms in the Luni object.
        """
        return self._file_loader.indices

    # @property
    # def names(self) -> NDArray[np.str_]:
    #     """Retrieve the names of the atoms in the Luni object.
    #
    #     Returns:
    #         NDArray[np.str_]: A NumPy array containing the names of the atoms in the Luni object.
    #     """
    #     return self.arc.atoms.names

    @property
    def names(self) -> NDArray[np.str_]:
        """Retrieve the names of the atoms in the Luni object.

        Returns:
            NDArray[np.str_]: A NumPy array containing the names of the atoms in the Luni object.
        """
        return self._file_loader.names

    # @property
    # def elements(self) -> NDArray[np.str_]:
    #     """Retrieve the elements of the atoms in the Luni object.
    #
    #     Returns:
    #         NDArray[np.str_]: A NumPy array containing the elements of the atoms in the Luni object.
    #     """
    #     return self.arc.atoms.elements

    @property
    def elements(self) -> NDArray[np.str_]:
        """Retrieve the elements of the atoms in the Luni object.

        Returns:
            NDArray[np.str_]: A NumPy array containing the elements of the atoms in the Luni object.
        """
        return self._file_loader.elements

    # @property
    # def types(self) -> NDArray[np.str_]:
    #     """Retrieve the types of the atoms in the Luni object.
    #
    #     Returns:
    #         NDArray[np.str_]: A NumPy array containing the types of the atoms in the Luni object.
    #     """
    #     return self.arc.atoms.types

    @property
    def coordinates(self) -> NDArray[np.float32]:
        """Retrieve the coordinates of the atoms in the Luni object.

        Returns:
            NDArray[np.float32]: A NumPy array containing the coordinates of the atoms in the Luni object.
        """
        return self._file_loader.coordinates

    @property
    def resnames(self) -> NDArray[np.str_]:
        """Retrieve the residue names of the atoms in the Luni object.

        Returns:
            NDArray[np.str_]: A NumPy array containing the residue names of the atoms in the Luni object.
        """
        return self._file_loader.resnames

    @property
    def resids(self) -> NDArray[np.int32]:
        """Retrieve the residue IDs of the atoms in the Luni object.

        Returns:
            NDArray[np.int32]: A NumPy array containing the residue IDs of the atoms in the Luni object.
        """
        return self._file_loader.resids

    # @property
    # def chainids(self) -> NDArray[np.int32]:
    #     """Retrieve the chain IDs of the atoms in the Luni object.
    #
    #     Returns:
    #         NDArray[np.str_]: A NumPy array containing the chain IDs of the atoms in the Luni object.
    #     """
    #     return self.arc.chains.ids

    @property
    def resindices(self):
        return self._file_loader.resindices

    # @property
    # def chainlabels(self) -> NDArray[np.str_]:
    #     """Retrieve the chain labels of the atoms in the Luni object.
    #
    #     Returns:
    #         NDArray[np.str_]: A NumPy array containing the chain labels of the atoms in the Luni object.
    #     """
    #     return self._file_loader.chainids

    # @property
    # def chainauths(self) -> NDArray[np.str_]:
    #     """Retrieve the chain auths of the atoms in the Luni object.
    #
    #     Returns:
    #         NDArray[np.str_]: A NumPy array containing the chain auths of the atoms in the Luni object.
    #     """
    #     return self._file_loader.auths

    @property
    def n_atoms(self) -> int:
        """Retrieve the number of atoms in the Luni object.

        Returns:
            int: The number of atoms in the Luni object.
        """
        return self._file_loader.n_atoms

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
