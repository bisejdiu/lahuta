"""
Placeholder for the universe module.
"""


import numpy as np
from MDAnalysis.core.topology import Topology

from lahuta.core.atom_assigner import AtomTypeAssigner
from lahuta.core.base import FileLoader
from lahuta.core.groups import AtomGroup
from lahuta.core.loaders import CIFLoader, PDBLoader
from lahuta.utils.atom_types import assign_radii, find_hydrogen_bonded_atoms


class Universe:
    def __init__(self, file_name=None, *args):
        if isinstance(file_name, Topology):
            raise NotImplementedError(
                "Initializing Universe from a Topology object is not supported."
            )

        file_loader = self._create_file_loader(file_name if file_name else args[0])
        self._universe = file_loader.load(*args)
        self.mol = file_loader.mol

        self._universe.atoms = AtomGroup(self._universe.atoms)
        self.atoms._u = self

        # self.hbond_array = HydrogenBondFinder(molecule_io.mol).find_hydrogen_bonded_atoms()
        self.hbond_array = find_hydrogen_bonded_atoms(self.mol)

        top_attr = {
            "resname": self.atoms.resnames,
            "name": self.atoms.names,
        }
        self.topology_attributes = top_attr
        atomtype_assigner = AtomTypeAssigner(
            self.mol, self.atoms, top_attr, legacy=True, parallel=False
        )
        atypes_array = atomtype_assigner.assign_atom_types()

        # save atypes_array to pickle file
        # import pickle

        # with open(
        #     "/home/bisejdiu/tutorials/lahuta-notebooks/data/res_atypes_array.pkl", "wb"
        # ) as f:
        #     pickle.dump(atypes_array, f)
        # self._extend_topology("vdw_radii", RadiusAssigner(molecule_io.mol).assign_radii())
        self._extend_topology("vdw_radii", assign_radii(self.mol))
        self._extend_topology("atom_types", atypes_array)

    @property
    def universe(self):
        return self

    @staticmethod
    def _create_file_loader(file_name: str) -> FileLoader:
        file_ext = file_name.split(".")[-1]
        if file_ext.lower() == "cif":
            return CIFLoader(file_name)
        else:
            return PDBLoader(file_name)

    def _extend_topology(self, attrname: str, values: np.ndarray):
        self.add_TopologyAttr(attrname, values)

    def select_atoms(self, *args, **kwargs) -> AtomGroup:
        return self.atoms.select_atoms(*args, **kwargs)

    def compute_neighbors(self, *args, **kwargs):
        return self.atoms.compute_neighbors(*args, **kwargs)

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


# class Universe(mda.Universe):
#     def __init__(self, topology=None, *args):

#         self.filename = None
#         if isinstance(topology, Topology):
#             raise NotImplementedError("Initializing Universe from a Topology object is not supported.")

#         file_loader = self.create_file_loader(topology)
#         self._universe = file_loader.load(*args)
#         self.mol = file_loader.mol

#         super().__init__(self._universe._topology)

#         self.atoms = AtomGroup(self.atoms)

#         # self.hbond_array = HydrogenBondFinder(molecule_io.mol).find_hydrogen_bonded_atoms()
#         self.hbond_array = find_hydrogen_bonded_atoms(self.mol)

#         top_attr = {
#             "resname": self.atoms.resnames,
#             "name": self.atoms.names,
#         }
#         atomtype_assigner = AtomTypeAssigner(self.mol, self.atoms, top_attr, legacy=False, parallel=False)
#         atypes_array = atomtype_assigner.assign_atom_types()

#         # self._extend_topology("vdw_radii", RadiusAssigner(molecule_io.mol).assign_radii())
#         self._extend_topology("vdw_radii", assign_radii(self.mol))
#         self._extend_topology("atom_types", atypes_array)

#     @staticmethod
#     def create_file_loader(file_name: str) -> FileLoader:
#         file_ext = file_name.split(".")[-1]
#         if file_ext.lower() == "cif":
#             return CIFLoader(file_name)
#         else: # file_ext.lower() == "pdb":
#             return PDBLoader(file_name)
#         # else:
#         #     raise ValueError(f"Unsupported file format: {file_ext}")


#     def _extend_topology(self, attrname: str, values: np.ndarray):
#         self.add_TopologyAttr(attrname, values)

#     def select_atoms(self, *args, **kwargs) -> AtomGroup:
#         return self.atoms.select_atoms(*args, **kwargs)

#     def compute_neighbors(self, *args, **kwargs):
#         return self.atoms.compute_neighbors(*args, **kwargs)

#     def __getattr__(self, attr):
#         # Delegate attribute access to the created universe
#         return getattr(self._universe, attr)

#     def __repr__(self):
#         return f"<Lahuta Universe with {self.atoms.n_atoms} atoms>"

#     def __str__(self):
#         return self.__repr__()
