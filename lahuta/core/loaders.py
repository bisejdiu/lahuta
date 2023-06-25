import itertools
from io import StringIO

import gemmi
import MDAnalysis as mda
import numpy as np
import pandas as pd
from openbabel import openbabel as ob

from lahuta.core.base import FileLoader
from lahuta.core.cra import Atoms, Chains, Residues
from lahuta.core.obmol import OBMol


class CIFLoader:
    def __init__(self, file_path):
        self.file_path = file_path
        self.block = gemmi.cif.read_file(file_path).sole_block()
        self.structure = gemmi.make_structure_from_block(self.block)
        self.atom_site_data = self.block.get_mmcif_category("_atom_site.")

        self.n_atoms = len(self.atom_site_data.get("Cartn_x"))

        self.chains = Chains(self.atom_site_data)
        self.residues = Residues(self.atom_site_data)
        self.atoms = Atoms(self.atom_site_data)

        self.atoms_df = pd.DataFrame(
            {
                "atom_name": self.atoms.names,
                "chain_name": self.chains.auths,
                "res_id": self.residues.resids,
                "res_name": self.residues.resnames,
            }
        )

        self.coords_array = self.extract_positions(self.atom_site_data)

    def extract_positions(self, atom_site_data):
        coords_array = np.zeros((self.n_atoms, 3))

        coords_array[:, 0] = atom_site_data.get("Cartn_x")
        coords_array[:, 1] = atom_site_data.get("Cartn_y")
        coords_array[:, 2] = atom_site_data.get("Cartn_z")

        return coords_array

    def _get_atom_index(self, atom_name, chain_name, res_id, res_name):
        index = self.atoms_df[
            (self.atoms_df["atom_name"] == atom_name)
            & (self.atoms_df["chain_name"] == chain_name)
            & (self.atoms_df["res_id"] == res_id)
            & (self.atoms_df["res_name"] == res_name)
        ].index

        if len(index) != 1:
            raise ValueError("Atom is not unique or does not exist")

        return int(index[0])

    def create_obmol(self):
        obmol = OBMol()

        added_residues = set()
        ob_res = None
        for idx, (chain, residue, atom) in enumerate(
            zip(self.chains, self.residues, self.atoms)
        ):
            _, chain_id = chain
            resname, resnumber, _ = residue
            _, atom_id, element = atom

            if ob_res is None or ob_res not in added_residues:
                ob_res = obmol.create_residue_obmol(resnumber, resname, chain_id)
                added_residues.add(ob_res)

            obmol.create_atom_obmol(
                int(atom_id), element, self.coords_array[idx], ob_res
            )

        obmol.perceive_bonds()
        for connection in self.structure.connections:
            prt1, prt2 = connection.partner1, connection.partner2
            atom1 = self._get_atom_index(
                prt1.atom_name, prt1.chain_name, prt1.res_id.seqid.num, prt1.res_id.name
            )
            atom2 = self._get_atom_index(
                prt2.atom_name, prt2.chain_name, prt2.res_id.seqid.num, prt2.res_id.name
            )

            obmol.create_bond_obmol(atom1, atom2)

        obmol.end_modify(True)

        return obmol.obmol

    def create_mda_universe(self):
        print(self.chains.mapping)
        resnames, resids, chain_ids = [], [], []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    resids.append(residue.seqid.num)
                    resnames.append(residue.name)
                    chain_ids.append(self.chains.mapping[chain.name])

        n_residues = len(resids)

        mda_universe = mda.Universe.empty(
            n_atoms=self.n_atoms,
            n_residues=n_residues,
            atom_resindex=self.residues.resindices,
            residue_segindex=chain_ids,
            trajectory=True,
        )

        mda_universe.add_TopologyAttr("names", self.atoms.names)
        mda_universe.add_TopologyAttr("type", self.atoms.elements)
        mda_universe.add_TopologyAttr("elements", self.atoms.elements)
        mda_universe.add_TopologyAttr("resnames", resnames)
        mda_universe.add_TopologyAttr("resids", resids)
        mda_universe.add_TopologyAttr("segids", np.array(["PROT"], dtype=object))

        mda_universe.atoms.positions = self.coords_array  # type: ignore

        return mda_universe

    def load(self):
        obmol = self.create_obmol()
        universe = self.create_mda_universe()
        return obmol, universe


class CIFLoader_old(FileLoader):
    def load(self, *args):
        self._load_obabel()
        # universe = self._create_mda_universe_from_obmol(self.mol)
        universe = mda.Universe(
            mda.lib.util.NamedStream(StringIO(self.pdb_str), "dummy.pdb")
        )
        return self.mol, universe

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
