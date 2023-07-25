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
        _file_loader = TopologyLoader.from_mda(args[0])  # type: ignore
        _mdag = _file_loader.to("mda")
        return _file_loader, _mdag

    def _initialize_from_files(self, files: str) -> Tuple[BaseLoader, AtomGroupType]:
        _file_loader = self._get_file_loader(files)
        _mdag = _file_loader.to("mda")
        return _file_loader, _mdag

    def _get_file_loader(self, *files: str) -> BaseLoader:
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
        # print("value size", values.size, values.shape)

        self._topattr_handler.init_topattr(attrname, attrname)
        self._mdag.universe.add_TopologyAttr(attrname, values)

    # def select_atoms(self, *args, **kwargs) -> mda.AtomGroup:
    #     return self.atoms.select_atoms(*args, **kwargs)

    def _build_atom_mapping(self, ag: AtomGroupType) -> NDArray[np.int64]:
        max_index = np.max(ag.universe.atoms.indices)
        atom_mapping = np.full(max_index + 1, -1, dtype=np.int64)
        atom_mapping[ag.indices] = np.arange(ag.n_atoms)
        return atom_mapping

    def ready(self) -> None:
        """
        Prepare instance for computations by transforming the molecule and assigning atom types.
        """

        assert self._file_loader is not None
        self._mol = self._file_loader.to("mol")
        self._mapping = self._build_atom_mapping(self.to("mda").universe.atoms)

        # TODO: remove array from the variable names by instead using type hints
        # self.hbond_array = find_hydrogen_bonded_atoms(self._mdag, self._mol)
        # print("...", hbond_array)
        # print("self._mdag", self._mdag, type(self._mdag))
        atomtype_assigner = AtomTypeAssigner(self._mdag, self._mol, self._mapping)
        ag_types = atomtype_assigner.assign_atom_types()
        og_atoms = self._mdag.universe.atoms

        # ag_types = AtomTypeAssigner(self._mol, self._mdag).assign_atom_types()
        reference_array = np.zeros((og_atoms.n_atoms, ag_types.shape[1]))
        # reference_hbond_array = np.zeros((og_atoms.n_atoms, 6), dtype=int)

        ix = self._mdag.indices
        full_ag_atypes = reference_array.copy()
        # full_ag_hbonds = reference_array[:, :6].copy()
        full_ag_atypes[ix] = ag_types
        # reference_hbond_array[ix] = hbond_array
        # self.hbond_array = reference_hbond_array

        # print("...2", self.hbond_array)

        # self._extend_topology("vdw_radii", v_radii_assignment(self.atoms.elements))
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

        Args:
        ----
        radius (float, optional): The cutoff radius. Default is 5.0.
        skip_adjacent (bool, optional): Whether to skip adjacent. Default is True.
        res_dif (int, optional): The residue difference to consider. Default is 1.

        Returns:
        -------
        NeighborPairs : np.ndarray
            An array of shape (n_atoms, n_neighbors) where each row contains the indices of the neighbors of the atom in the row.

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
        """Retrieve the file format from a file name.

        Args:
            file_name (str): The name of the file.

        Returns:
            tuple: The file format and a boolean indicating if it is pdb or pdb.gz.
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

        Args:
        ----
        fmt (str): The format to convert to. Currently supported formats are "mda" and "mol".
        args (list): Trajectory file name(s). Only required if converting to "mda".

        Returns:
        -------
        Universe : Universe
            A new Universe instance in the specified format.
        """
        if fmt == "mol" and self._mol is not None:
            return self._mol
        # if fmt == "mda":
        #     return self._mdag
        return self._file_loader.to(fmt)

    @property
    def arc(self) -> Union[None, ARC]:
        return self._file_loader.arc

    def __repr__(self) -> str:
        return f"<Lahuta Universe with {self._mdag.n_atoms} atoms>"

    def __str__(self) -> str:
        return self.__repr__()
