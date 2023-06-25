import tempfile
from abc import ABC, abstractmethod
from typing import Optional

import gemmi
from openbabel import openbabel as ob


class ModelProcessor:
    def __init__(self, mol=None):
        self.chain_numbs = {}
        self.atom_model_map = {}
        self.mol = ob.OBMol() if mol is None else mol
        # self.mol.SetChainsPerceived()
        # self.mol.BeginModify()  # type: ignore
        # self.mol.SetChainsPerceived()

    def process_structure(self, structure):
        for model in structure:
            self.process_model(model)
        return self.mol, self.atom_model_map

    def process_model(self, model):
        for chain in model:
            self.process_chain(chain, model)

    def process_chain(self, chain, model):
        chain_name = chain.name
        if chain_name not in self.chain_numbs:
            self.chain_numbs[chain_name] = len(self.chain_numbs) + 1

        for residue in chain:
            self.process_residue(residue, chain_name, model)

    def process_residue(self, residue, chain_name, model):
        ob_res = self._create_residue_obmol(
            self.mol, residue, self.chain_numbs[chain_name]
        )

        for atom in residue:
            self.process_atom(atom, residue, chain_name, model, ob_res)

    def process_atom(self, atom, residue, chain_name, model, ob_res):
        atom_address = gemmi.AtomAddress(
            chain_name, residue.seqid, residue.name, atom.name
        )
        self.atom_model_map[str(atom_address)] = model

        ob_atom = self._create_atom_obmol(self.mol, atom)
        self._add_atoms_to_residue(ob_res, ob_atom)

    def _create_residue_obmol(self, mol, residue, chain_name):
        res_id = residue.seqid.num
        res_name = residue.name

        ob_res = mol.NewResidue()
        ob_res.SetChainNum(chain_name)
        ob_res.SetNum(str(res_id))
        ob_res.SetName(str(res_name))

        return ob_res

    def _create_atom_obmol(self, mol, atom):
        ob_atom = mol.NewAtom(atom.serial - 1)
        ob_atom.SetAtomicNum(gemmi.Element(atom.element.name).atomic_number)
        ob_atom.SetVector(atom.pos.x, atom.pos.y, atom.pos.z)
        ob_atom.SetFormalCharge(atom.charge)

        return ob_atom

    def _add_atoms_to_residue(self, ob_res, ob_atom):
        ob_res.AddAtom(ob_atom)
        ob_res.SetHetAtom(
            ob_atom, not gemmi.find_tabulated_residue(ob_res.GetName()).is_standard()
        )
        ob_res.SetSerialNum(ob_atom, ob_res.GetSerialNum(ob_atom))

    @staticmethod
    def create_obmol_bond(mol, atom1_id, atom2_id):
        ob_atom1 = mol.GetAtomById(atom1_id)
        ob_atom2 = mol.GetAtomById(atom2_id)

        for neighbor in ob.OBAtomAtomIter(ob_atom1):
            if neighbor.GetId() == ob_atom2.GetId():
                continue

        ob_bond = mol.NewBond()
        ob_bond.SetBegin(ob_atom1)
        ob_bond.SetEnd(ob_atom2)
        ob_bond.SetBondOrder(1)  # TODO: Research how to set bond order

    @staticmethod
    def perceive_obmol_properties(mol):
        mol.SetChainsPerceived()
        # mol.SetAromaticPerceived()
        # mol.SetAtomTypesPerceived()
        # mol.SetChiralityPerceived()
        # mol.SetRingTypesPerceived()
        # mol.SetPartialChargesPerceived()

        return mol


class FileLoader(ABC):
    """
    Base class for file loaders.

    Parameters
    ----------
    file_name : str
        The name of the file to load.

    in_format : str, optional
        The format of the input file. If not specified, the format will be
        inferred from the file extension.
    """

    def __init__(self, file_name: str, in_format: Optional[str] = None):
        self.file_name = file_name
        self.in_format = (
            in_format if in_format is not None else self.file_name.split(".")[-1]
        )
        self.mol = None
        self.pdb_str = None

    @abstractmethod
    def load(self, *args):
        """Load the file and return the universe"""
        raise NotImplementedError("Subclasses must implement this method")

    def read_cif_as_pdb(self, cif_file: str):
        """Convert a CIF file to a PDB file. Uses MDAnalysis to do the conversion."""

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pdb", delete=False
        ) as pdb_file:
            log_level = ob.cvar.obErrorLog.GetOutputLevel()
            ob.cvar.obErrorLog.SetOutputLevel(0)
            ob_conv = ob.OBConversion()
            ob_conv.SetInFormat("cif")
            temp_mol = ob.OBMol()
            ob_conv.ReadFile(temp_mol, cif_file)

            ob_conv.SetOutFormat("pdb")
            ob_conv.WriteFile(temp_mol, pdb_file.name)
            pdb_str = ob_conv.WriteString(temp_mol)

            pdb_file.write(pdb_str)
            pdb_file.flush()

            ob_conv = ob.OBConversion()
            ob_conv.SetInFormat("pdb")
            mol = ob.OBMol()
            ob_conv.ReadString(mol, pdb_str)

            ob.cvar.obErrorLog.SetOutputLevel(log_level)

            return pdb_str, mol

    def _load_obabel(self):
        if self.in_format.lower() == "cif":
            pdb_str, self.mol = self.read_cif_as_pdb(self.file_name)
            self.pdb_str = pdb_str
        else:
            ob_conv = ob.OBConversion()
            ob_conv.SetInFormat(self.in_format)
            mol = ob.OBMol()
            ob_conv.ReadFile(mol, self.file_name)
            self.mol = mol
