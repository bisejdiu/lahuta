import gemmi
import pandas as pd
from openbabel import openbabel as ob  # type: ignore


class OBMol:
    def __init__(self):
        self.mol = None

    def create_residue_obmol(self, resid, resname, chain_id):
        assert self.mol is not None, "Molecule is not initialized"
        ob_res = self.mol.NewResidue()
        ob_res.SetChainNum(int(chain_id))
        ob_res.SetNum(str(resid))
        ob_res.SetName(str(resname))
        return ob_res

    def create_atom_obmol(
        self, idx, atom_name, atom_id, atom_element, atom_pos, ob_residue
    ):
        assert self.mol is not None, "Molecule is not initialized"
        ob_atom = self.mol.NewAtom(idx)
        ob_atom.SetAtomicNum(gemmi.Element(atom_element).atomic_number)
        ob_atom.SetType(atom_name)
        ob_atom.SetVector(float(atom_pos[0]), float(atom_pos[1]), float(atom_pos[2]))
        ob_atom.SetFormalCharge(0)  # atom.charge)

        self.add_atoms_to_residue(ob_atom, ob_residue)

    def create_bond_obmol(self, atom1_id, atom2_id):
        assert self.mol is not None, "Molecule is not initialized"
        ob_atom1 = self.mol.GetAtomById(atom1_id)
        ob_atom2 = self.mol.GetAtomById(atom2_id)

        for neighbor in ob.OBAtomAtomIter(ob_atom1):
            if neighbor.GetId() == ob_atom2.GetId():
                continue

        ob_bond = self.mol.NewBond()
        ob_bond.SetBegin(ob_atom1)
        ob_bond.SetEnd(ob_atom2)
        # ob_bond.SetBondOrder(1)  # TODO: Research how to set bond order

    def add_atoms_to_residue(self, ob_atom, ob_res):
        ob_res.AddAtom(ob_atom)
        ob_res.SetHetAtom(
            ob_atom, not gemmi.find_tabulated_residue(ob_res.GetName()).is_standard()
        )
        ob_res.SetSerialNum(ob_atom, ob_res.GetSerialNum(ob_atom))

    def perceive_bonds(self):
        if self.mol:
            self.mol.ConnectTheDots()
            self.mol.PerceiveBondOrders()

    def perceive_properties(self):
        # self.mol.SetChainsPerceived()
        if self.mol:
            self.mol.SetAromaticPerceived()
            self.mol.SetAtomTypesPerceived()
            self.mol.SetChiralityPerceived()
            self.mol.SetRingTypesPerceived()
            self.mol.SetPartialChargesPerceived()

        return self.mol

    def end_modify(self, nuke_perceived_data: bool = True):
        if self.mol:
            self.mol.EndModify(nuke_perceived_data)

    def create_mol(self, arc, connections=None):
        chains, residues, atoms = arc.chains, arc.residues, arc.atoms
        coords = atoms.coordinates
        if connections is None:
            connections = []

        atoms_df = pd.DataFrame(
            {
                "atom_name": arc.atoms.names,
                "chain_name": arc.chains.auths,
                "res_id": arc.residues.resids,
                "res_name": arc.residues.resnames,
            }
        )

        self.mol = ob.OBMol()
        self.mol.BeginModify()

        ob_res = None
        added_residues = set()
        for idx, (chain, residue, atom) in enumerate(zip(chains, residues, atoms)):
            atom_name, atom_id, element = atom["name"], atom["id"], atom["element"]
            resname, resid = residue["resname"], residue["resid"]
            chain_id = chain["id"]

            _cra_ = (chain_id, resid, resname)
            if ob_res is None or _cra_ not in added_residues:
                ob_res = self.create_residue_obmol(resid, resname, chain_id)
                added_residues.add(_cra_)

            self.create_atom_obmol(
                idx, atom_name, int(atom_id), element, coords[idx], ob_res
            )

        self.perceive_bonds()
        for connection in connections:
            prt1, prt2 = connection.partner1, connection.partner2
            atom1 = self._get_atom_index(
                atoms_df,
                prt1.atom_name,
                prt1.chain_name,
                prt1.res_id.seqid.num,
                prt1.res_id.name,
            )
            atom2 = self._get_atom_index(
                atoms_df,
                prt2.atom_name,
                prt2.chain_name,
                prt2.res_id.seqid.num,
                prt2.res_id.name,
            )

            self.create_bond_obmol(atom1, atom2)

        self.perceive_properties()

        self.end_modify(True)
        self.mol.SetChainsPerceived()  # type: ignore

    def _get_atom_index(self, atoms_df, atom_name, chain_name, res_id, res_name):
        index = atoms_df[
            (atoms_df["atom_name"] == atom_name)
            & (atoms_df["chain_name"] == chain_name)
            & (atoms_df["res_id"] == res_id)
            & (atoms_df["res_name"] == res_name)
        ].index

        if len(index) != 1:
            raise ValueError("Atom is not unique or does not exist")

        return int(index[0])
