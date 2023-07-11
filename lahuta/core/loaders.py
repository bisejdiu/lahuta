import gemmi
import MDAnalysis as mda
import numpy as np
import pandas as pd

from lahuta.core.base import FileLoader
from lahuta.core.cra import Atoms, Chains, Residues
from lahuta.core.obmol import OBMol


class GemmiLoader:
    def __init__(self, file_path, is_pdb=False):
        self.file_path = file_path
        if is_pdb:
            self.structure = gemmi.read_pdb(self.file_path)
            self.block = self.structure.make_mmcif_document().sole_block()
        else:
            self.block = gemmi.cif.read(self.file_path).sole_block()
            self.structure = gemmi.make_structure_from_block(self.block)

        self.atom_site_data = self.block.get_mmcif_category("_atom_site.")

        self._chains = Chains().from_gemmi(self.atom_site_data)
        self._residues = Residues().from_gemmi(self.atom_site_data)
        self._atoms = Atoms().from_gemmi(self.atom_site_data)

    @property
    def n_atoms(self):
        return self._atoms.names.shape[0]

    @property
    def chains(self):
        return self._chains

    @property
    def residues(self):
        return self._residues

    @property
    def atoms(self):
        return self._atoms


class MDALoader:
    def __init__(self, file_path, traj_path=None):
        self.file_path = file_path
        self.universe = mda.Universe(self.file_path)

        self._chains = Chains().from_mda(self.universe)
        self._residues = Residues().from_mda(self.universe)
        self._atoms = Atoms().from_mda(self.universe)

    @property
    def n_atoms(self):
        return self._atoms.names.shape[0]

    @property
    def chains(self):
        return self._chains

    @property
    def residues(self):
        return self._residues

    @property
    def atoms(self):
        return self._atoms


class StructureLoader:
    def __init__(self, file_path, is_pdb=False):
        self.file_path = file_path
        # extension = file_path.split(".")[-1]
        self.loader = GemmiLoader(file_path, is_pdb=is_pdb)

        self.atoms_df = pd.DataFrame(
            {
                "atom_name": self.loader.atoms.names,
                "chain_name": self.loader.chains.auths,
                "res_id": self.loader.residues.resids,
                "res_name": self.loader.residues.resnames,
            }
        )

        self.coords_array = self.extract_positions(self.loader.atom_site_data)

    def extract_positions(self, atom_site_data):
        coords_array = np.zeros((self.loader.n_atoms, 3))

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

        ob_res = None
        added_residues = set()
        for idx, (chain, residue, atom) in enumerate(
            zip(self.loader.chains, self.loader.residues, self.loader.atoms)
        ):
            _, chain_id = chain
            resname, resnumber, _ = residue
            atom_name, atom_id, element = atom

            cra = (chain_id, resnumber, resname)
            if ob_res is None or cra not in added_residues:
                # print(cra)
                ob_res = obmol.create_residue_obmol(resnumber, resname, chain_id)
                added_residues.add(cra)

            obmol.create_atom_obmol(
                atom_name, int(atom_id), element, self.coords_array[idx], ob_res
            )

        obmol.perceive_bonds()
        for connection in self.loader.structure.connections:
            prt1, prt2 = connection.partner1, connection.partner2
            atom1 = self._get_atom_index(
                prt1.atom_name, prt1.chain_name, prt1.res_id.seqid.num, prt1.res_id.name
            )
            atom2 = self._get_atom_index(
                prt2.atom_name, prt2.chain_name, prt2.res_id.seqid.num, prt2.res_id.name
            )

            obmol.create_bond_obmol(atom1, atom2)

        obmol.perceive_properties()

        obmol.end_modify(True)
        obmol.mol.SetChainsPerceived()  # type: ignore

        return obmol.mol

    def create_mda_universe(self):
        # TODO: Add icodes and ids to the universe
        resnames, resids, chain_ids = [], [], []
        for model in self.loader.structure:
            for chain in model:
                for residue in chain:
                    resids.append(residue.seqid.num)
                    resnames.append(residue.name)
                    chain_ids.append(self.loader.chains.mapping[chain.name])

        n_residues = len(resids)

        mda_universe = mda.Universe.empty(
            n_atoms=self.loader.n_atoms,
            n_residues=n_residues,
            atom_resindex=self.loader.residues.resindices,
            residue_segindex=chain_ids,
            trajectory=True,
        )

        mda_universe.add_TopologyAttr("names", self.loader.atoms.names)
        mda_universe.add_TopologyAttr("type", self.loader.atoms.types)
        mda_universe.add_TopologyAttr("elements", self.loader.atoms.elements)
        mda_universe.add_TopologyAttr("resnames", resnames)
        mda_universe.add_TopologyAttr("resids", resids)
        mda_universe.add_TopologyAttr("segids", np.array(["PROT"], dtype=object))

        mda_universe.atoms.positions = self.coords_array  # type: ignore

        return mda_universe

    def load(self):
        obmol = self.create_obmol()
        universe = self.create_mda_universe()
        return obmol, universe


class PDBLoader(FileLoader):
    def load(self, *args):
        self._load_obabel()
        universe = mda.Universe(self.file_name, *args)
        return self.mol, universe
