"""Generate an Open Babel molecule from an ARC object.

Generating an Open Babel molecule from an ARC object is done by creating a new OBMol object and calling the
create_mol method with the ARC object as an argument. This is necessary because of the quite involved 
process of creating a molecule using the provided Open Babel API. It, unfortunately, does not support 
vectorized operations, so we have to create the molecule atom by atom. There is a performance penalty
for this. 

Classes:
    ```
    OBMol: A class for generating an Open Babel molecule from an ARC object.
    ```

Example:
    ``` py
    obmol = OBMol()
    obmol.create_mol(arc)
    obmol.mol
    ```

"""
from typing import Any, Optional

import gemmi
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from openbabel import openbabel as ob

from lahuta._types.openbabel import MolAtomType, MolResType, MolType, MolTypeWrapper

from . import ARC


class OBMol:
    """A class for generating an Open Babel molecule from an ARC object."""

    def __init__(self) -> None:
        self.mol: Optional[MolType] = None

    def create_residue_obmol(
        self,
        resid: NDArray[np.int32],
        resname: NDArray[np.str_],
        chain_id: NDArray[np.int32],
    ) -> MolResType:
        """Create a new residue in the molecule.

        Args:
            resid (NDArray[np.int32]): The residue ID.
            resname (NDArray[np.str_]): The residue name.
            chain_id (NDArray[np.int32]): The chain ID.

        Returns:
            MolResType: A new residue instance.
        """
        assert self.mol is not None, "Molecule is not initialized"
        ob_res = self.mol.NewResidue()
        ob_res.SetChainNum(int(chain_id))
        ob_res.SetNum(str(resid))
        ob_res.SetName(str(resname))
        return ob_res

    def create_atom_obmol(
        self,
        idx: int,
        atom_name: str,
        atom_element: str,
        atom_pos: NDArray[np.float32],
        ob_residue: MolResType,
    ) -> None:
        """Create a new atom in the molecule.

        Args:
            idx (int): The index of the atom.
            atom_name (str): The name of the atom.
            atom_element (str): The type of atom to add.
            atom_pos (NDArray[np.float_]): The position of the atom.
            ob_residue (MolResType): The residue the atom belongs to.
        """
        assert self.mol is not None, "Molecule is not initialized"
        ob_atom: MolAtomType = self.mol.NewAtom(idx)
        ob_atom.SetAtomicNum(gemmi.Element(atom_element).atomic_number)
        ob_atom.SetType(atom_name)
        ob_atom.SetVector(float(atom_pos[0]), float(atom_pos[1]), float(atom_pos[2]))
        ob_atom.SetFormalCharge(0)  # atom.charge)

        self.add_atoms_to_residue(ob_atom, ob_residue)

    def create_bond_obmol(self, atom1_id: int, atom2_id: int) -> None:
        """Create a bond between two atoms in the molecule.

        Args:
            atom1_id (int): The ID of the first atom.
            atom2_id (int): The ID of the second atom.
        """
        assert self.mol is not None, "Molecule is not initialized"
        ob_atom1 = self.mol.GetAtomById(atom1_id)
        ob_atom2 = self.mol.GetAtomById(atom2_id)

        for neighbor in ob.OBAtomAtomIter(ob_atom1):
            if neighbor.GetId() == ob_atom2.GetId():
                continue

        ob_bond = self.mol.NewBond()
        ob_bond.SetBegin(ob_atom1)
        ob_bond.SetEnd(ob_atom2)

    def add_atoms_to_residue(self, ob_atom: MolAtomType, ob_res: MolResType) -> None:
        """Add an atom to a residue.

        Args:
            ob_atom (MolAtomType): The atom to add.
            ob_res (MolResType): The residue to add the atom to.
        """
        ob_res.AddAtom(ob_atom)
        ob_res.SetHetAtom(ob_atom, not gemmi.find_tabulated_residue(ob_res.GetName()).is_standard())  # type: ignore
        ob_res.SetSerialNum(ob_atom, ob_res.GetSerialNum(ob_atom))

    def perceive_bonds(self) -> None:
        """Identify all the bonds in the molecule."""
        if self.mol:
            self.mol.ConnectTheDots()
            log_level = ob.cvar.obErrorLog.GetOutputLevel()
            ob.cvar.obErrorLog.SetOutputLevel(0)

            self.mol.PerceiveBondOrders()
            ob.cvar.obErrorLog.SetOutputLevel(log_level)

    def perceive_properties(self) -> MolType | None:
        """Identify properties of the molecule.

        Returns:
            MolType | None: The molecule with its properties perceived.
        """
        if self.mol:
            self.mol.SetAromaticPerceived()
            self.mol.SetAtomTypesPerceived()
            self.mol.SetChiralityPerceived()
            self.mol.SetRingTypesPerceived()
            self.mol.SetPartialChargesPerceived()

        return self.mol

    def end_modify(self, nuke_perceived_data: bool = True) -> None:
        """End the modification of the molecule.

        Args:
            nuke_perceived_data (bool, optional): If True, perceived data is deleted. Default is True.
        """
        if self.mol:
            self.mol.EndModify(nuke_perceived_data)

    def create_mol(self, arc: ARC, connections: Optional[Any] = None) -> None:  # noqa: ANN401
        """Create a new molecule from an ARC (Atomic Record Collection) object.

        Args:
            arc (ARC): An object representing atomic data.
            connections (Optional[Any], optional): A list of connections between atoms. Default is None.
        """
        chains, residues, atoms = arc.chains, arc.residues, arc.atoms
        coords = atoms.coordinates

        atoms_df = pd.DataFrame(
            {
                "atom_name": arc.atoms.names,
                "chain_name": arc.chains.auths,
                "res_id": arc.residues.resids,
                "res_name": arc.residues.resnames,
            }
        )

        # use MolTypeWrapper to create a new molecule
        self.mol = MolTypeWrapper(ob.OBMol()).mol
        self.mol.BeginModify()

        ob_res = None
        added_residues: set[tuple[NDArray[Any], NDArray[Any], NDArray[Any]]] = set()
        for idx, (chain, residue, atom) in enumerate(zip(chains, residues, atoms, strict=True)):
            atom_id = int(atom["id"])
            atom_name, element = atom["name"], atom["element"]
            resname: NDArray[np.str_] = residue["resname"]
            resid: NDArray[np.int32] = residue["resid"]
            chain_id: NDArray[np.int32] = chain["id"]

            _cra_ = (chain_id, resid, resname)
            if ob_res is None or _cra_ not in added_residues:
                ob_res = self.create_residue_obmol(resid, resname, chain_id)
                added_residues.add(_cra_)

            self.create_atom_obmol(atom_id, str(atom_name), str(element), coords[idx], ob_res)

        self.perceive_bonds()

        if connections is None:
            connections = []
        assert connections is not None
        for conn_data in connections:
            atom1: int = self._get_atom_index(
                atoms_df,
                conn_data.partner1.atom_name,
                conn_data.partner1.chain_name,
                conn_data.partner1.res_id_seq_num,
                conn_data.partner1.res_id_name,
            )

            atom2: int = self._get_atom_index(
                atoms_df,
                conn_data.partner2.atom_name,
                conn_data.partner2.chain_name,
                conn_data.partner2.res_id_seq_num,
                conn_data.partner2.res_id_name,
            )

            self.create_bond_obmol(atom1, atom2)

        self.perceive_properties()

        self.end_modify(True)
        self.mol.SetChainsPerceived()

    def _get_atom_index(
        self,
        atoms_df: pd.DataFrame,
        atom_name: str,
        chain_name: str,
        res_id: int,
        res_name: str,
    ) -> int:
        """Get the index of an atom in a DataFrame.

        Args:
            atoms_df (pd.DataFrame): A DataFrame with atomic data.
            atom_name (str): The name of the atom.
            chain_name (str): The name of the chain.
            res_id (int): The residue ID.
            res_name (str): The name of the residue.

        Returns:
            int: The index of the atom in the DataFrame.
        """
        index = atoms_df[
            (atoms_df["atom_name"] == atom_name)
            & (atoms_df["chain_name"] == chain_name)
            & (atoms_df["res_id"] == res_id)
            & (atoms_df["res_name"] == res_name)
        ].index

        if not len(index):
            raise ValueError("Atom is not unique or does not exist")

        return int(index[0])
