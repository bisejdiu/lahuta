import gemmi
from openbabel import openbabel as ob


class OBMol:
    def __init__(self):
        self.obmol = ob.OBMol()
        self.obmol.SetChainsPerceived()
        self.obmol.BeginModify()

    def create_residue_obmol(self, resid, resname, chain_id):
        ob_res = self.obmol.NewResidue()
        ob_res.SetChainNum(int(chain_id))
        ob_res.SetNum(str(resid))
        ob_res.SetName(str(resname))
        return ob_res

    def create_atom_obmol(self, atom_id, atom_element, atom_pos, ob_residue):
        ob_atom = self.obmol.NewAtom(atom_id - 1)
        ob_atom.SetAtomicNum(gemmi.Element(atom_element).atomic_number)
        ob_atom.SetVector(float(atom_pos[0]), float(atom_pos[1]), float(atom_pos[2]))
        # ob_atom.SetFormalCharge(atom.charge)

        self.add_atoms_to_residue(ob_atom, ob_residue)

    def create_bond_obmol(self, atom1_id, atom2_id):
        ob_atom1 = self.obmol.GetAtomById(atom1_id)
        ob_atom2 = self.obmol.GetAtomById(atom2_id)

        for neighbor in ob.OBAtomAtomIter(ob_atom1):
            if neighbor.GetId() == ob_atom2.GetId():
                continue

        ob_bond = self.obmol.NewBond()
        ob_bond.SetBegin(ob_atom1)
        ob_bond.SetEnd(ob_atom2)
        ob_bond.SetBondOrder(1)  # TODO: Research how to set bond order

    def add_atoms_to_residue(self, ob_atom, ob_res):
        ob_res.AddAtom(ob_atom)
        ob_res.SetHetAtom(
            ob_atom, not gemmi.find_tabulated_residue(ob_res.GetName()).is_standard()
        )
        ob_res.SetSerialNum(ob_atom, ob_res.GetSerialNum(ob_atom))

    def perceive_bonds(self):
        self.obmol.ConnectTheDots()
        self.obmol.PerceiveBondOrders()

    def perceive_properties(self):
        self.obmol.SetChainsPerceived()
        self.obmol.SetAromaticPerceived()
        self.obmol.SetAtomTypesPerceived()
        self.obmol.SetChiralityPerceived()
        self.obmol.SetRingTypesPerceived()
        self.obmol.SetPartialChargesPerceived()

        return self.obmol

    def end_modify(self, nuke_perceived_data: bool = True):
        self.obmol.EndModify(nuke_perceived_data)
