"""
Placeholder for the universe module.
"""


from typing import Any

import MDAnalysis as mda
import numpy as np
from MDAnalysis.core.universe import (_check_file_like,
                                      _generate_from_topology,
                                      _resolve_coordinates,
                                      _topology_from_file_like)
from openbabel import openbabel as ob

from lahuta.core.atom_assigner import AtomTypeAssigner
from lahuta.core.groups import AtomGroup
from lahuta.utils.atom_types import assign_radii, find_hydrogen_bonded_atoms


class MoleculeIO:
    def __init__(self, file_name: str, coordinates: Any = (), in_format: str = None):
        self.coordinates = coordinates
        self.file_name = file_name
        self.in_format = in_format
        self.mol = self._create_obmol()

    def _create_obmol(self):
        if self.in_format is None:
            try:
                self.in_format = self.file_name.split(".")[-1]
            except AttributeError as exc:
                raise AttributeError("You must either specify the input format using the in_format argument or pass a file name with a file extension.") from exc

        log_level = ob.cvar.obErrorLog.GetOutputLevel()
        ob.cvar.obErrorLog.SetOutputLevel(0)

        ob_conv = ob.OBConversion()
        ob_conv.SetInFormat(self.in_format)
        mol = ob.OBMol()
        ob_conv.ReadFile(mol, self.file_name)
        ob.cvar.obErrorLog.SetOutputLevel(log_level)

        return mol

    def create_mda_universe_from_obmol(self):
        """
        Create an MDAnalysis Universe from an Open Babel molecule object.
        """
        n_atoms = self.mol.NumAtoms()
        n_residues = self.mol.NumResidues()
        
        # Initialize topology attributes in a dictionary and coordinates array
        topology_attributes = {
            'name': np.empty(n_atoms, dtype=object),
            'type': np.empty(n_atoms, dtype=object),
            'resname': np.empty(n_residues, dtype=object),
            'resid': np.empty(n_residues, dtype=int),
            'element': np.empty(n_atoms, dtype=object),
            'segid': np.array(['PROT'], dtype=object)
        }
        coords = np.empty((n_atoms, 3))
        resindices = np.empty(n_atoms, dtype=int)

        atom_counter = 0
        for res_idx, res in enumerate(ob.OBResidueIter(self.mol)):
            for atom in ob.OBResidueAtomIter(res):
                topology_attributes['name'][atom_counter] = res.GetAtomID(atom)
                topology_attributes['type'][atom_counter] = atom.GetType()
                topology_attributes['element'][atom_counter] = ob.GetSymbol(atom.GetAtomicNum())
                coords[atom_counter] = [atom.GetX(), atom.GetY(), atom.GetZ()]
                resindices[atom_counter] = res_idx
                atom_counter += 1
            
            topology_attributes['resname'][res_idx] = res.GetName()
            topology_attributes['resid'][res_idx] = res.GetIdx() + 1  # 1-based indexing

        # Create the Universe
        universe = mda.Universe.empty(n_atoms=n_atoms,
                                    n_residues=n_residues,
                                    atom_resindex=resindices,
                                    residue_segindex=np.array([0] * n_residues, dtype=int),
                                    trajectory=True)

        # Add topology attributes
        for attr, values in topology_attributes.items():
            # print ('Adding topology attribute:', attr, 'with values:', values, 'to universe.')
            universe.add_TopologyAttr(attr, values)

        # Add positions
        if self.coordinates:
            universe.load_new(self.coordinates)
        else:
            universe.atoms.positions = coords

        return universe, topology_attributes


class Universe(mda.Universe):
    def __init__(self, topology=None, *args, **kwargs):

        self.filename = None

        coordinates = _resolve_coordinates(self.filename, *args, format=None, all_coordinates=False)

        molecule_io = MoleculeIO(topology, coordinates, **kwargs)
        universe, topology_attributes = molecule_io.create_mda_universe_from_obmol()
        self.mol = molecule_io.mol

        super().__init__(universe._topology)

        self.topology_attributes = topology_attributes
        self.trajectory = universe.trajectory
        self.atoms.positions = universe.atoms.positions
        self.atoms = AtomGroup(self.atoms)
        # self.hbond_array = HydrogenBondFinder(molecule_io.mol).find_hydrogen_bonded_atoms()
        self.hbond_array = find_hydrogen_bonded_atoms(self.mol)

        atomtype_assigner = AtomTypeAssigner(molecule_io.mol, self.atoms, self.topology_attributes, legacy=False, parallel=False)
        atypes_array = atomtype_assigner.assign_atom_types()

        # self._extend_topology("vdw_radii", RadiusAssigner(molecule_io.mol).assign_radii())
        self._extend_topology("vdw_radii", assign_radii(self.mol))
        self._extend_topology("atom_types", atypes_array)

    def _extend_topology(self, attrname: str, values: np.ndarray):
        self.add_TopologyAttr(attrname, values)

    def select_atoms(self, *args, **kwargs) -> AtomGroup:
        return self.atoms.select_atoms(*args, **kwargs)

    def compute_neighbors(self, *args, **kwargs):
        return self.atoms.compute_neighbors(*args, **kwargs)

    def __repr__(self):
        return f"<Lahuta Universe with {self.atoms.n_atoms} atoms>"

    def __str__(self):
        return self.__repr__()

