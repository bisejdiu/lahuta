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
from typing_extensions import Self

from lahuta._types.mdanalysis import AtomGroupType
from lahuta.config.atom_types import BASE_AA_CONVERSION, RESIDUE_SYNONYMS
from lahuta.config.defaults import GEMMI_SUPPRTED_FORMATS, MDA_SUPPORTED_FORMATS
from lahuta.core import loader as L
from lahuta.core.neighbors import NeighborPairs
from lahuta.core.topology import LahutaCPPLoader, TopologyLoader
from lahuta.core.topology.loaders import load_file
from lahuta.lib._lahuta import LahutaCPP
from lahuta.utils.array_utils import cross_interaction_indices

IntArray, StrArray, FloatArray, AnyArray = list[int], list[str], list[float], list[Any]
if TYPE_CHECKING:
    from numpy.typing import NDArray

    IntArray = NDArray[np.int32]
    StrArray = NDArray[np.str_]
    FloatArray = NDArray[np.float32]
    AnyArray = NDArray[Any]

    from lahuta._types.mdanalysis import TrajectoryType
    from lahuta._types.openbabel import MolType

__all__ = ["Luni"]


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
        # structure: "str | Path | AtomGroupType",
        structure: str,
        trajectories: str | list[str] | None = None,
        b_iso_name: str = "tempfactor",
    ) -> None:
        fmts: str | set[str] = ""
        self._data: LahutaCPP | AtomGroupType
        self._input_structure: str | Path | "AtomGroupType" = structure
        structure = str(structure) if isinstance(structure, Path) else structure
        # print("->", structure, trajectories)
        match (structure, trajectories):
            case (mda.AtomGroup(atoms=s), None):
                # print("converting from MDA")
                self._data = TopologyLoader.from_mda(s)  # type: ignore
            case (str(s), None):
                # Check if we can use GemmiLoader.
                file_format, is_pdb = Luni._check_gemmi_support(s)
                if file_format:
                    # print("-> to bench")
                    self._data = LahutaCPPLoader(s).luni
                    self._fd = L.LoaderFactory.load(s, file_format)
                    print("->", self._fd.to_ir())
                elif s.upper().split(".")[-1] in MDA_SUPPORTED_FORMATS:
                    # print("MDA")
                    self._data = TopologyLoader((s)).to("mda")
                else:
                    fmts = self._get_supported_fmts()
                    fmts = ", ".join(fmts)
                    raise ValueError(f"Unsupported format for structure: {s}! \nSupported formats are: {fmts}.")
            case (str(s), str(t)):
                # If trajectories are provided, use TopologyLoader
                self._data = TopologyLoader(s, t).to("mda")
            case _:
                fmts = self._get_supported_fmts()
                fmts = ", ".join(fmts)
                raise ValueError(
                    f"Invalid input! {structure=} and {trajectories=} are not valid inputs. \nSupported formats are: {fmts}."
                )

        self._mda = L.MDAnalysisLoader.from_ir(self._fd.to_ir())
        self._structure, self._trajectories = structure, trajectories

        self.b_iso_name = b_iso_name
        if self.b_iso_name != "tempfactor":
            self.extend_topology(self.b_iso_name, self._mda.universe.atoms.tempfactors)

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
        """Compute the neighbors of each atom in the Luni object."""

        if image is not None and backend != "gemmi":
            raise ValueError("Image filtering is only supported with the 'gemmi' backend!")

        mda = self._mda
        if target_spec is not None:
            union_indices = np.union1d(self.indices, target_spec.indices)
            mda = self._mda.universe.atoms[union_indices]

        neighbors = self._data.find_neighbors(radius, res_dif)
        pairs, distances = neighbors.get_pairs(), np.array(neighbors.get_distances_sq())

        if target_spec is not None:
            cross_indices = cross_interaction_indices(pairs, self.indices, target_spec.indices)
            pairs, distances = pairs[cross_indices], distances[cross_indices]

        ns = NeighborPairs(self)
        ns.set_neighbors(pairs, distances)
        ns._ns = neighbors
        return ns

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
        # if fmt == "mol":
        #     self._mol = self._data.to(fmt)
        return getattr(self, f"_{fmt}")  # type: ignore

    @property
    def indices(self) -> IntArray:
        """Retrieve the indices of the atoms in the Luni object.

        Returns:
            NDArray[np.int32]: A NumPy array containing the indices of the atoms in the Luni object.
        """
        return self._data.indices

    @property
    def ids(self) -> IntArray:
        """Retrieve the indices of the atoms in the Luni object.

        Returns:
            NDArray[np.int32]: A NumPy array containing the indices of the atoms in the Luni object.
        """
        return self._data.indices

    @property
    def names(self) -> StrArray:
        """Retrieve the names of the atoms in the Luni object.

        Returns:
            NDArray[np.str_]: A NumPy array containing the names of the atoms in the Luni object.
        """
        return self._data.names

    @property
    def elements(self) -> StrArray:
        """Retrieve the elements of the atoms in the Luni object.

        Returns:
            NDArray[np.str_]: A NumPy array containing the elements of the atoms in the Luni object.
        """
        return self._data.elements

    @property
    def coordinates(self) -> FloatArray:
        """Retrieve the coordinates of the atoms in the Luni object.

        Returns:
            NDArray[np.float32]: A NumPy array containing the coordinates of the atoms in the Luni object.
        """
        return self._data.positions

    @property
    def resnames(self) -> StrArray:
        """Retrieve the residue names of the atoms in the Luni object.

        Returns:
            NDArray[np.str_]: A NumPy array containing the residue names of the atoms in the Luni object.
        """
        return self._data.resnames

    @property
    def resids(self) -> IntArray:
        """Retrieve the residue IDs of the atoms in the Luni object.

        Returns:
            NDArray[np.int32]: A NumPy array containing the residue IDs of the atoms in the Luni object.
        """
        return self._data.resids

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
        return self._data.resindices

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
        return self._data.n_atoms

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

    def read_file(self, file_path: str) -> LahutaCPP | AtomGroupType:
        """Read the file using the provided loader.

        Args:
            file_path (str): The path to the file.
            loader (BaseLoader): The loader to use for reading the file.
        """
        return load_file(file_path)

    @staticmethod
    def _get_supported_fmts() -> set[str]:
        return GEMMI_SUPPRTED_FORMATS.union({x.lower() for x in MDA_SUPPORTED_FORMATS})

    def extend_topology(self, attrname: str, values: AnyArray) -> None:
        """Add new topology attributes to the Luni.

        Args:
            attrname (str): The name of the attribute.
            values (NDArray[Any]): The values of the attribute.
        """
        from .topology.topattrs import AtomAttrClassHandler

        topattr_handler = AtomAttrClassHandler()
        topattr_handler.init_topattr(attrname, attrname)
        self._mda.universe.add_TopologyAttr(attrname, values)
