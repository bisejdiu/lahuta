"""
Placeholder for the universe module.
"""

from typing import Literal

import MDAnalysis as mda
import numpy as np
from MDAnalysis.core.topology import Topology

from lahuta.config.defaults import GEMMI_SUPPRTED_FORMATS
from lahuta.config.smarts import AVAILABLE_ATOM_TYPES
from lahuta.core._loaders import GemmiLoader, TopologyLoader
from lahuta.core.atom_assigner import AtomTypeAssigner
from lahuta.core.groups import AtomGroup
from lahuta.core.neighbor_finder import NeighborSearch
from lahuta.core.topattrs import AtomAttrClassHandler
from lahuta.utils.atom_types import find_hydrogen_bonded_atoms, v_radii_assignment


# class Universe:
#     def __init__(self, file_name=None, *args):
#         if isinstance(file_name, Topology):
#             raise NotImplementedError(
#                 "Initializing Universe from a Topology object is not supported."
#             )
class Universe:
    def __init__(self, *args):
        if len(args) == 0:
            raise ValueError(
                "Must pass either a file name or an MDAnalysis Universe or MDAnalysis.AtomGroup instance"
            )

        if isinstance(args[0], mda.Universe) or isinstance(args[0], mda.AtomGroup):
            if len(args) > 1:
                raise ValueError(
                    "When passing an MDAnalysis.Universe or MDAnalysis.AtomGroup instance, no other arguments are allowed"
                )
            self._initialize_from_universe(args[0])
        else:
            self._initialize_from_files(args)

        # self.atoms._u = self

        self.atoms = self.file_loader.atoms
        self.residues = self.file_loader.residues
        self.chains = self.file_loader.chains
        self.universe = self.file_loader.to("mda")

        self.mol = None
        self.hbond_array = None
        self._ready = False
        self._topattr_handler = AtomAttrClassHandler()

        # print("number of atoms: ", self.atoms.n_atoms)

    def _initialize_from_universe(self, uniatom):
        self.file_loader = TopologyLoader.from_mda(
            uniatom
        )  # use self._file_loader to get the correct loader
        self._universe = self.file_loader.to("mda")
        self._universe.atoms = AtomGroup(self._universe.atoms[uniatom.indices])
        # self.atoms = self._universe.atoms[uniatom.indices]
        # self._universe.atoms = self._universe.atoms[uniatom.indices]
        # self._universe._topology = uniatom.universe._topology
        # self.file_loader = TopologyLoader.from_mda(uniatom)
        # self._initialize_common()

    def _initialize_from_files(self, files):
        for file in files:
            if not isinstance(file, str):
                raise ValueError(
                    "All arguments must be filenames when not providing an MDAnalysis Universe or MDAnalysis.AtomGroup instance"
                )
        # Assuming the first file is the topology file
        self.file_loader = self._file_loader(files[0])
        # self._universe = self.file_loader.to("mda", *files)
        self.universe = self.file_loader.to("mda", *files)

    @classmethod
    def from_mda(cls, mda_universe):
        return cls(mda_universe)

    # @property
    # def universe(self):
    #     return self

    @staticmethod
    def _file_loader(file_name: str):  # -> FileLoader:
        file_format, is_pdb = Universe.get_format(file_name)
        if file_format is not None:
            return GemmiLoader(file_name, is_pdb=is_pdb)

        return TopologyLoader(file_name)

    def _extend_topology(self, attrname: str, values: np.ndarray):
        self._topattr_handler.init_topattr(attrname, attrname)
        self.universe.add_TopologyAttr(attrname, values)

    # def select_atoms(self, *args, **kwargs) -> mda.AtomGroup:
    #     return self.atoms.select_atoms(*args, **kwargs)

    def ready(self):
        """
        Prepare instance for computations by transforming the molecule and assigning atom types.
        """

        self.mol = self.file_loader.to("mol")

        # TODO: remove array from the variable names by instead using type hints
        self.hbond_array = find_hydrogen_bonded_atoms(self.mol)
        atomtype_assigner = AtomTypeAssigner(self.mol, self.universe.atoms)
        atypes_array = atomtype_assigner.assign_atom_types()

        self._extend_topology("vdw_radii", v_radii_assignment(self.atoms.elements))
        for atom_type in AVAILABLE_ATOM_TYPES:
            self._extend_topology(
                atom_type.name.lower(), atypes_array[:, atom_type.value]
            )

        self._ready = True

    # TODO: rename to find_neighbors
    def compute_neighbors(
        self,
        radius=5.0,
        ignore_hydrogens=True,
        res_dif=1,
    ):
        """
        Compute the neighbors of each atom in the Universe.

        Args:
        ----
        radius (float, optional): The cutoff radius. Default is 5.0.
        ignore_hydrogens (bool, optional): Whether to ignore hydrogens. Default is True.
        skip_adjacent (bool, optional): Whether to skip adjacent. Default is True.
        res_dif (int, optional): The residue difference to consider. Default is 1.

        Returns:
        -------
        NeighborPairs : np.ndarray
            An array of shape (n_atoms, n_neighbors) where each row contains the indices of the neighbors of the atom in the row.

        """
        if not self._ready:
            self.ready()

        neighbors = NeighborSearch(self)
        return neighbors.compute(
            radius=radius,
            ignore_hydrogens=ignore_hydrogens,
            res_dif=res_dif,
        )

    @staticmethod
    def get_format(file_name):
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

    def to(self, fmt: Literal["mda", "mol"], *args):
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
        return self.file_loader.to(fmt, *args)

    # def __getattr__(self, attr):
    #     # Delegate attribute access to the created universe
    #     return getattr(self._universe, attr)

    # def __dir__(self):
    #     universe_dir = set(dir(self._universe))
    #     self_dir = set(super().__dir__())
    #     return sorted(self_dir.union(universe_dir))

    def __repr__(self):
        return f"<Lahuta Universe with {self.universe.atoms.n_atoms} atoms>"

    def __str__(self):
        return self.__repr__()
