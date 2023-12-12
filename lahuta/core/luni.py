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
from lahuta.core.assigners import AtomTypeAssigner
from lahuta.core.neighbors import BaseNeighborSearch, NeighborPairs
from lahuta.core.neighbors.backends import GemmiNeighborSearch, MDAnalysisNeighborSearch
from lahuta.core.topology import GemmiLoader, TopologyLoader
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
        self, structure: Union[str, Path, "AtomGroupType"], trajectories: Optional[str | list[str]] = None
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
        from .topology.topattrs import AtomAttrClassHandler

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

        # neighbors = MDAnalysisNeighborSearch(self.to("mda"))
        neighbors: BaseNeighborSearch
        if backend == "gemmi":
            import gemmi
            # TODO(bisejdiu): image is not being passed
            structure = gemmi.read_pdb(self._input_structure)
            neighbors = GemmiNeighborSearch(mda, structure)
        elif backend == "mda":
            neighbors = MDAnalysisNeighborSearch(mda)
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
        """Filter atoms from this Luni using a selection string to a new Luni.

        `lahuta` uses the MDAnalysis selection language to filter atoms. It defines
        several additional keywords in addition to the standard MDAnalysis keywords.

        The documentation below is taken and lightly adabpted from the MDAnalysis
        documentation. See https://userguide.mdanalysis.org/stable/selections.html
        for more.

        The selection parser understands the following CASE SENSITIVE *keywords*:

        **Simple selections**

            protein, backbone, nucleic, water, ion, lipid, etc.
                selects all atoms that belong to a standard set of residues;
                a protein is identfied by a hard-coded set of residue names so
                it  may not work for esoteric residues. TODO (bisejdiu): add
                a list of all supported keywords belonging to this category.
            segid *seg-name*
                select by segid (as given in the topology), e.g. ``segid 4AKE``
                or ``segid DMPC``
            resid *residue-number-range*
                resid can take a single residue number or a range of numbers. A
                range consists of two numbers separated by a colon (inclusive)
                such as ``resid 1:5``. A residue number ("resid") is taken
                directly from the topology.
                If icodes are present in the topology, then these will be
                taken into account.  Ie 'resid 163B' will only select resid
                163 with icode B while 'resid 163' will select only residue 163.
                Range selections will also respect icodes, so 'resid 162-163B'
                will select all residues in 162 and those in 163 up to icode B.
            resnum *resnum-number-range*
                resnum is the canonical residue number; typically it is set to
                the residue id in the original PDB structure.
            resname *residue-name*
                select by residue name, e.g. ``resname LYS``
            name *atom-name*
                select by atom name (as given in the topology). Often, this is
                force field dependent. Example: ``name CA`` (for C&alpha; atoms)
                or ``name OW`` (for SPC water oxygen)
            type *atom-type*
                select by atom type; this is either a string or a number and
                depends on the force field; it is read from the topology file
                (e.g. the CHARMM PSF file contains numeric atom types). It has
                non-sensical values when a PDB or GRO file is used as a topology
            atom *seg-name*  *residue-number*  *atom-name*
                a selector for a single atom consisting of segid resid atomname,
                e.g. ``DMPC 1 C2`` selects the C2 carbon of the first residue of
                the DMPC segment
            altloc *alternative-location*
                a selection for atoms where alternative locations are available,
                which is often the case with high-resolution crystal structures
                e.g. `resid 4 and resname ALA and altloc B` selects only the
                atoms of ALA-4 that have an altloc B record.
            moltype *molecule-type*
                select by molecule type, e.g. ``moltype Protein_A``. At the
                moment, only the TPR format defines the molecule type.
            record_type *record_type*
                for selecting either ATOM or HETATM from PDB-like files.
                e.g. ``select_atoms('name CA and not record_type HETATM')``
            smarts *SMARTS-query*
                select atoms using Daylight's SMARTS queries, e.g. ``smarts
                [#7;R]`` to find nitrogen atoms in rings. Requires RDKit.
                All matches (max 1000) are combined as a unique match
            chiral *R | S*
                select a particular stereocenter. e.g. ``name C and chirality
                S`` to select only S-chiral carbon atoms.  Only ``R`` and
                ``S`` will be possible options but other values will not raise
                an error.

        **Boolean**

            not
                all atoms not in the selection, e.g. ``not protein`` selects
                all atoms that aren't part of a protein
            and, or
                combine two selections according to the rules of boolean
                algebra, e.g. ``protein and not resname ALA LYS``
                selects all atoms that belong to a protein, but are not in a
                lysine or alanine residue

        **Geometric**

            around *distance*  *selection*
                selects all atoms a certain cutoff away from another selection,
                e.g. ``around 3.5 protein`` selects all atoms not belonging to
                protein that are within 3.5 Angstroms from the protein
            point *x* *y* *z*  *distance*
                selects all atoms within a cutoff of a point in space, make sure
                coordinate is separated by spaces,
                e.g. ``point 5.0 5.0 5.0  3.5`` selects all atoms within 3.5
                Angstroms of the coordinate (5.0, 5.0, 5.0)
            prop [abs] *property*  *operator*  *value*
                selects atoms based on position, using *property*  **x**, **y**,
                or **z** coordinate. Supports the **abs** keyword (for absolute
                value) and the following *operators*: **<, >, <=, >=, ==, !=**.
                For example, ``prop z >= 5.0`` selects all atoms with z
                coordinate greater than 5.0; ``prop abs z <= 5.0`` selects all
                atoms within -5.0 <= z <= 5.0.
            sphzone *radius* *selection*
                Selects all atoms that are within *radius* of the center of
                geometry of *selection*
            sphlayer *inner radius* *outer radius* *selection*
                Similar to sphzone, but also excludes atoms that are within
                *inner radius* of the selection COG
            cyzone *externalRadius* *zMax* *zMin* *selection*
                selects all atoms within a cylindric zone centered in the
                center of geometry (COG) of a given selection,
                e.g. ``cyzone 15 4 -8 protein and resid 42`` selects the
                center of geometry of protein and resid 42, and creates a
                cylinder of external radius 15 centered on the COG. In z, the
                cylinder extends from 4 above the COG to 8 below. Positive
                values for *zMin*, or negative ones for *zMax*, are allowed.
            cylayer *innerRadius* *externalRadius* *zMax* *zMin* *selection*
                selects all atoms within a cylindric layer centered in the
                center of geometry (COG) of a given selection,
                e.g. ``cylayer 5 10 10 -8 protein`` selects the center of
                geometry of protein, and creates a cylindrical layer of inner
                radius 5, external radius 10 centered on the COG. In z, the
                cylinder extends from 10 above the COG to 8 below. Positive
                values for *zMin*, or negative ones for *zMax*, are allowed.

        **Connectivity**

            byres *selection*
                selects all atoms that are in the same segment and residue as
                selection, e.g. specify the subselection after the byres keyword
            bonded *selection*
                selects all atoms that are bonded to selection
                eg: ``select name H and bonded name O`` selects only hydrogens
                bonded to oxygens

        **Index**

            bynum *index-range*
                selects all atoms within a range of (1-based) inclusive indices,
                e.g. ``bynum 1`` selects the first atom in the universe;
                ``bynum 5:10`` selects atoms 5 through 10 inclusive. All atoms
                in the :class:`~MDAnalysis.core.universe.Universe` are
                consecutively numbered, and the index runs from 1 up to the
                total number of atoms.
            index *index-range*
                selects all atoms within a range of (0-based) inclusive indices,
                e.g. ``index 0`` selectsssor.results the first atom in the universe;
                ``index 5:10`` selects atoms 6 through 11 inclusive. All atoms
                in the :class:`~MDAnalysis.core.universe.Universe` are
                consecutively numbered, and the index runs from 0 up to the
                total number of atoms - 1.

        **Preexisting selections**

            group `group-name`
                selects the atoms in the :class:`AtomGroup` passed to the
                function as a keyword argument named `group-name`. Only the
                atoms common to `group-name` and the instance
                :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms`
                was called from will be considered, unless ``group`` is
                preceded by the ``global`` keyword. `group-name` will be
                included in the parprocesing just by comparison of atom indices.
                This means that it is up to the user to make sure the
                `group-name` group was defined in an appropriate
                :class:`~MDAnalysis.core.universe.Universe`.
            global *selection*
                by default, when issuing
                :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` from an
                :class:`~MDAnalysis.core.groups.AtomGroup`, selections and
                subselections are returned intersected with the atoms of that
                instance. Prefixing a selection term with ``global`` causes its
                selection to be returned in its entirety.  As an example, the
                ``global`` keyword allows for
                ``lipids.select_atoms("around 10 global protein")`` --- where
                ``lipids`` is a group that does not contain any proteins. Were
                ``global`` absent, the result would be an empty selection since
                the ``protein`` subselection would itself be empty. When issuing
                :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` from a
                :class:`~MDAnalysis.core.universe.Universe`, ``global`` is
                ignored.

        **Dynamic selections**
            Not tested & not documented & probably not working.

        Args:
            selection (str): The selection string to use for filtering.

        Returns:
            Luni: A new Luni instance containing the filtered atoms.
        """
        return self.__class__(self._mda.select_atoms(selection))
        import pandas as pd

        mda_copy = self._mda.select_atoms(selection).copy()

        uv = mda.Universe.empty(
            n_atoms=mda_copy.n_atoms,
            n_residues=mda_copy.n_residues,
            n_segments=mda_copy.n_segments,
            atom_resindex=pd.factorize(mda_copy.resindices)[0],
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

        assert uv.atoms is not None
        uv.atoms.positions = mda_copy.positions
        uv.filename = mda_copy.universe.filename

        return self.__class__(uv.atoms)

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
