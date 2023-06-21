import itertools

import MDAnalysis as mda
import numpy as np
from openbabel import openbabel as ob

from lahuta.core.base import FileLoader


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
            resindices[atom_counter] = residue.GetNum()

            topology_attributes["resname"][res_idx] = residue.GetName()
            topology_attributes["resid"][
                res_idx
            ] = residue.GetNum()  # res_idx + 1  # 1-based indexing

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
        universe.atoms.positions = coords  # type: ignore

        return universe


class PDBLoader(FileLoader):
    def load(self, *args):
        self._load_obabel()
        universe = mda.Universe(self.file_name, *args)
        return universe


# class MoleculeIO:
#     def __init__(self, file_loader: FileLoader, coordinates: Any = ()):
#         self.coordinates = coordinates
#         self.file_loader = file_loader
#         self.mol, self.universe = self._create_obmol_and_mda_universe()

#     def _create_obmol_and_mda_universe(self):
#         mol, universe = self.file_loader.load()
#         if self.coordinates:
#             universe.load_new(self.coordinates)
#         return mol, universe
