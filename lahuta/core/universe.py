"""
Placeholder for the universe module.
"""

from pathlib import Path

import MDAnalysis as mda
import numpy as np
from openbabel import openbabel as ob

from lahuta.utils.cif_converter import convert_cif_to_pdb, read_cif_as_pdb

from ..utils.atom_types import (assign_atom_types, assign_radii,
                                find_hydrogen_bonded_atoms)
from .groups import AtomGroup
from .obabel import OBMol


def create_mda_universe_from_obmol(mol):
    n_atoms = mol.NumAtoms()
    n_residues = mol.NumResidues()
    
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
    for res_idx, res in enumerate(ob.OBResidueIter(mol)):
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
                                   trajectory=True)

    # Add topology attributes
    for attr, values in topology_attributes.items():
        print ('Adding topology attribute:', attr, 'with values:', values, 'to universe.')
        universe.add_TopologyAttr(attr, values)

    # Add positions
    universe.atoms.positions = coords

    return universe

class Universe(mda.Universe):
    """A subclass of MDAnalysis Universe that adds additional attributes and methods."""

    def __init__(self, *args, **kwargs):
        """A subclass of MDAnalysis Universe that adds some extra functionality."""
        self.mol = OBMol(*args)
        universe = create_mda_universe_from_obmol(self.mol)
        super().__init__(universe._topology, **kwargs)

        self.trajectory = universe.trajectory
        self.atoms.positions = universe.atoms.positions
        self.atoms = AtomGroup(self.atoms)
        self.hbond_array = find_hydrogen_bonded_atoms(self.mol)

        self._extend_topology("vdw_radii", assign_radii(self.mol))
        self._extend_topology("atom_types", assign_atom_types(self.mol, self.atoms))

    def _extend_topology(self, attrname: str, values: np.ndarray):
        """Extend the topology with a new attribute.

        Parameters
        ----------
        attrname : str
            The name of the attribute.

        values : np.ndarray
            The values of the attribute.
        """
        self.add_TopologyAttr(attrname, values)

    def select_atoms(self, *args, **kwargs) -> AtomGroup:
        """Select atoms.

        Wrapper around MDAnalysis Universe.select_atoms that returns a LahutaAtomGroup.
        """
        return self.atoms.select_atoms(*args, **kwargs)  # type: ignore

    def compute_neighbors(self, *args, **kwargs):
        """Compute neighbors.
        See Also
        --------
        :meth:`Lahuta.core.groups.AtomGroup.compute_neighbors`
        """
        return self.atoms.compute_neighbors(*args, **kwargs)

    def __repr__(self):
        return f"<Lahuta Universe with {self.atoms.n_atoms} atoms>"

    def __str__(self):
        return self.__repr__()

# class Universe(mda.Universe):
#     """A subclass of MDAnalysis Universe that adds additional attributes and methods."""

#     def __init__(self, *args, **kwargs):
#         """A subclass of MDAnalysis Universe that adds some extra functionality."""
#         # get the first argument from the args list
#         self.filename = args[0]

#         suffix = "pdb"
#         if isinstance(self.filename, str):
#             suffix = self.filename.split(".")[-1]
#         elif isinstance(self.filename, Path):
#             suffix = self.filename.suffix[1:]
#             self.filename = str(self.filename)

#         if suffix == "cif":
#             from io import StringIO

#             from lahuta.utils.cif_converter import (convert_cif_to_pdb,
#                                                     read_cif_as_pdb)

#             # pdb_str, mol = convert_cif_to_pdb(self.filename)
#             pdb_str, mol = read_cif_as_pdb(self.filename)
#             super().__init__(
#                 mda.lib.util.NamedStream(StringIO(pdb_str), "dummy.pdb"), **kwargs
#             )
#             self.mol = mol
#         elif suffix == "pdb":

#             super().__init__(*args, **kwargs)
#             self.mol = OBMol(str(self.filename))

#         self.atoms = AtomGroup(self.atoms)
#         self.hbond_array = find_hydrogen_bonded_atoms(self.mol)

#         self._extend_topology("vdw_radii", assign_radii(self.mol))
#         self._extend_topology("atom_types", assign_atom_types(self.mol, self.atoms))

#     def _extend_topology(self, attrname: str, values: np.ndarray):
#         """Extend the topology with a new attribute.

#         Parameters
#         ----------
#         attrname : str
#             The name of the attribute.

#         values : np.ndarray
#             The values of the attribute.
#         """
#         self.add_TopologyAttr(attrname, values)

#     def select_atoms(self, *args, **kwargs) -> AtomGroup:
#         """Select atoms.

#         Wrapper around MDAnalysis Universe.select_atoms that returns a LahutaAtomGroup.
#         """
#         return self.atoms.select_atoms(*args, **kwargs)  # type: ignore

#     def compute_neighbors(self, *args, **kwargs):
#         """Compute neighbors.
#         See Also
#         --------
#         :meth:`Lahuta.core.groups.AtomGroup.compute_neighbors`
#         """
#         return self.atoms.compute_neighbors(*args, **kwargs)

#     def __repr__(self):
#         return f"<Lahuta Universe with {self.atoms.n_atoms} atoms>"

#     def __str__(self):
#         return self.__repr__()
