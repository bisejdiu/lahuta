"""
Placeholder for the universe module.
"""


import MDAnalysis as mda
import numpy as np
from MDAnalysis.core.topology import Topology
from MDAnalysis.lib.nsgrid import FastNS

from lahuta.config.defaults import GEMMI_SUPPRTED_FORMATS
from lahuta.config.smarts import AVAILABLE_ATOM_TYPES
from lahuta.core.atom_assigner import AtomTypeAssigner
from lahuta.core.base import FileLoader
from lahuta.core.loaders import CIFLoader, PDBLoader
from lahuta.core.neighbors import NeighborPairs
from lahuta.core.topattrs import AtomAttrClassHandler
from lahuta.utils.atom_types import assign_radii, find_hydrogen_bonded_atoms
from lahuta.utils.mda import mda_psuedobox_from_atomgroup


class Universe:
    def __init__(self, file_name=None, *args):
        if isinstance(file_name, Topology):
            raise NotImplementedError(
                "Initializing Universe from a Topology object is not supported."
            )

        file_loader = self._create_file_loader(file_name if file_name else args[0])
        self.mol, self._universe = file_loader.load()

        self.atoms._u = self

        self.hbond_array = find_hydrogen_bonded_atoms(self.mol)

        top_attr = {
            "resname": self.atoms.resnames,
            "name": self.atoms.names,
        }
        self.topology_attributes = top_attr
        atomtype_assigner = AtomTypeAssigner(
            self.mol, self.atoms, top_attr, legacy=False, parallel=False
        )
        atypes_array = atomtype_assigner.assign_atom_types()

        self._topattr_handler = AtomAttrClassHandler()

        self._extend_topology("vdw_radii", assign_radii(self.mol))
        # self._extend_topology("atom_types", atypes_array)
        for atom_type in AVAILABLE_ATOM_TYPES:
            self._extend_topology(
                atom_type.name.lower(), atypes_array[:, atom_type.value]
            )

    @property
    def universe(self):
        return self

    @staticmethod
    def _create_file_loader(file_name: str):  # -> FileLoader:
        file_format, is_pdb = Universe.get_format(file_name)
        if file_format is not None:
            return CIFLoader(file_name, is_pdb=is_pdb)
        else:
            # TODO: Channel to an MDA loader
            print("Not Supported Format")
            return PDBLoader(file_name)

    def _extend_topology(self, attrname: str, values: np.ndarray):
        self._topattr_handler.init_topattr(attrname, attrname)
        self.add_TopologyAttr(attrname, values)

    # def select_atoms(self, *args, **kwargs) -> mda.AtomGroup:
    #     return self.atoms.select_atoms(*args, **kwargs)

    # def compute_neighbors(self, *args, **kwargs):
    #     return self.atoms.compute_neighbors(*args, **kwargs)

    def compute_neighbors(
        self,
        radius=5.0,
        ignore_hydrogens=True,
        # hbonds=False,
        skip_adjacent=True,
        res_dif=1,
    ):
        """Compute the neighbors of each atom in the Universe.

        Parameters
        ----------
        radius : float, optional
            The cutoff radius. Default is 5.0.
        ignore_hydrogens : bool, optional
            Whether to ignore hydrogens. Default is True.
        hbonds : bool, optional
            Whether to include hydrogen bonds. Default is False.

        Returns
        -------
        neighbors : np.ndarray
            An array of shape (n_atoms, n_neighbors) where each row contains the
            indices of the neighbors of the atom in the row.
        """
        if ignore_hydrogens:
            atomgroup = self.select_atoms("not name H*")
        else:
            atomgroup = self.atoms

        pairs, distances = self.get_neighbors(atomgroup, radius)

        if skip_adjacent:
            idx = self._remove_adjacent_residue_pairs(pairs, res_dif=res_dif)
            pairs = pairs[idx]
            distances = distances[idx]

        return NeighborPairs(self.atoms, pairs, distances)

    def get_neighbors(self, atomgroup=None, radius=5.0):
        """Get the neighbors of an atomgroup.

        Parameters
        ----------
        atomgroup : AtomGroup
            The atomgroup to get the neighbors of.

        Returns
        -------
        pairs : np.ndarray
            An array of shape (n_pairs, 2) where each row contains the indices of
            the atoms in the pair.
        distances : np.ndarray
            An array of shape (n_pairs,) containing the distances of each pair.
        """
        if atomgroup is None:
            atomgroup = self.atoms

        shift_coords, pseudobox = mda_psuedobox_from_atomgroup(atomgroup)

        gridsearch = FastNS(
            cutoff=radius, coords=shift_coords, box=pseudobox, pbc=False
        )
        neighbors = gridsearch.self_search()

        return atomgroup[neighbors.get_pairs()].indices, neighbors.get_pair_distances()

    def _remove_adjacent_residue_pairs(self, pairs, res_dif=1):
        """Remove pairs where the difference in residue ids is less than `res_dif`.

        Parameters
        ----------
        pairs : np.ndarray
            An array of shape (n_pairs, 2) where each row is a pair of atom indices.
        res_dif : int, optional
            The difference in residue ids to remove. Default is 1.

        Returns
        -------
        pairs : np.ndarray
            An array of shape (n_pairs, 2) where each row is a pair of atom indices.

        """
        resids = self.atoms.resids[pairs]
        return np.any(np.abs(resids - resids[:, ::-1]) > res_dif, axis=1)

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

    def __getattr__(self, attr):
        # Delegate attribute access to the created universe
        return getattr(self._universe, attr)

    def __dir__(self):
        universe_dir = set(dir(self._universe))
        self_dir = set(super().__dir__())
        return sorted(self_dir.union(universe_dir))

    def __repr__(self):
        return f"<Lahuta Universe with {self.atoms.n_atoms} atoms>"

    def __str__(self):
        return self.__repr__()
