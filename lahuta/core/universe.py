"""
Placeholder for the universe module.
"""

import itertools
from typing import Any

import MDAnalysis as mda
import numpy as np
from MDAnalysis.core.topology import Topology
from openbabel import openbabel as ob

from lahuta.core.atom_assigner import AtomTypeAssigner
from lahuta.core.base import FileLoader
from lahuta.core.groups import AtomGroup
from lahuta.utils.atom_types import assign_radii, find_hydrogen_bonded_atoms


class CIFLoader(FileLoader):
    def load(self, *args):
        self._load_obabel()
        universe = self._create_mda_universe_from_obmol(self.mol)
        return universe

    @staticmethod
    def _create_mda_universe_from_obmol(mol):
        """
        Create an MDAnalysis Universe from an Open Babel molecule object.
        """
        n_atoms = mol.NumAtoms()
        n_residues = mol.NumResidues()

        # Initialize topology attributes in a dictionary and coordinates array
        topology_attributes = {
            "name": np.empty(n_atoms, dtype=object),
            "type": np.empty(n_atoms, dtype=object),
            "resname": np.empty(n_residues, dtype=object),
            "resid": np.empty(n_residues, dtype=int),
            "element": np.empty(n_atoms, dtype=object),
            "segid": np.array(["PROT"], dtype=object),
        }
        coords = np.empty((n_atoms, 3))
        resindices = np.empty(n_atoms, dtype=int)

        atom_iterator = itertools.chain.from_iterable(
            ob.OBResidueAtomIter(res) for res in ob.OBResidueIter(mol)
        )
        for atom_counter, atom in enumerate(atom_iterator):
            residue = atom.GetResidue()
            res_idx = residue.GetIdx()
            element = ob.GetSymbol(atom.GetAtomicNum())
            topology_attributes["name"][atom_counter] = residue.GetAtomID(atom)
            topology_attributes["type"][atom_counter] = element
            topology_attributes["element"][atom_counter] = element
            coords[atom_counter] = [atom.GetX(), atom.GetY(), atom.GetZ()]
            resindices[atom_counter] = res_idx

            topology_attributes["resname"][res_idx] = residue.GetName()
            topology_attributes["resid"][res_idx] = res_idx + 1  # 1-based indexing

        # atom_counter = 0
        # for res_idx, res in enumerate(ob.OBResidueIter(mol)):
        #     for atom in ob.OBResidueAtomIter(res):
        #         element = ob.GetSymbol(atom.GetAtomicNum())
        #         topology_attributes["name"][atom_counter] = res.GetAtomID(atom)
        #         topology_attributes["type"][atom_counter] = element
        #         topology_attributes["element"][atom_counter] = element
        #         coords[atom_counter] = [atom.GetX(), atom.GetY(), atom.GetZ()]
        #         resindices[atom_counter] = res_idx
        #         atom_counter += 1

        # topology_attributes["resname"][res_idx] = res.GetName()
        # topology_attributes["resid"][res_idx] = res.GetIdx() + 1  # 1-based indexing

        # Create the Universe
        universe = mda.Universe.empty(
            n_atoms=n_atoms,
            n_residues=n_residues,
            atom_resindex=resindices,
            residue_segindex=np.array([0] * n_residues, dtype=int),
            trajectory=True,
        )

        # Add topology attributes
        for attr, values in topology_attributes.items():
            # print ('Adding topology attribute:', attr, 'with values:', values, 'to universe.')
            universe.add_TopologyAttr(attr, values)

        # Add positions
        universe.atoms.positions = coords

        return universe


class PDBLoader(FileLoader):
    def load(self, *args):
        self._load_obabel()
        universe = mda.Universe(self.file_name, *args)
        return universe


class MoleculeIO:
    def __init__(self, file_loader: FileLoader, coordinates: Any = ()):
        self.coordinates = coordinates
        self.file_loader = file_loader
        self.mol, self.universe = self._create_obmol_and_mda_universe()

    def _create_obmol_and_mda_universe(self):
        mol, universe = self.file_loader.load_molecule()
        if self.coordinates:
            universe.load_new(self.coordinates)
        return mol, universe


class Universe:
    def __init__(self, topology=None, *args):
        if isinstance(topology, Topology):
            raise NotImplementedError(
                "Initializing Universe from a Topology object is not supported."
            )

        file_loader = self._create_file_loader(topology)
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
