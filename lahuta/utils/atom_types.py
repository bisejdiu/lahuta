"""
Placehoder for the atom types and radii.
"""
import numpy as np

from openbabel import openbabel as ob

from ..config.atoms import PROT_ATOM_TYPES
from ..config.smarts import ATOM_TYPES
from ..config.residues import STANDARD_RESIDUES


def assign_atom_types(mol, atomgroup):
    """
    Assign atom types to each atom in the molecule.
    Atom types are defined in `ATOM_TYPES`
    """
    atypes = {
        "hbond acceptor": 0,
        "pos ionisable": 1,
        "carbonyl oxygen": 2,
        "weak hbond donor": 3,
        "carbonyl carbon": 4,
        "weak hbond acceptor": 5,
        "hbond donor": 6,
        "neg ionisable": 7,
        "aromatic": 8,
        "xbond acceptor": 9,
        "hydrophobe": 10,
    }

    atypes_array = np.zeros((mol.NumAtoms(), len(atypes)))
    for atom_type, smartsdict in ATOM_TYPES.items():

        for smarts in smartsdict.values():

            ob_smart = ob.OBSmartsPattern()
            ob_smart.Init(str(smarts))
            ob_smart.Match(mol)

            matches = [x[0] for x in ob_smart.GetMapList()]
            for match in matches:
                atom = mol.GetAtom(match)

                atypes_array[atom.GetId(), atypes[atom_type]] = 1

    # ALL WATER MOLECULES ARE HYDROGEN BOND DONORS AND ACCEPTORS
    for atom in atomgroup.select_atoms(
        "resname SOL HOH TIP3 TIP4 WAT W and not name H*"
    ):
        atypes_array[atom.index, atypes["hbond acceptor"]] = 1
        atypes_array[atom.index, atypes["hbond donor"]] = 1

    # OVERRIDE PROTEIN ATOM TYPING FROM DICTIONARY
    for residue in atomgroup.select_atoms(
        "resname " + " ".join(STANDARD_RESIDUES)
    ).residues:
        for atom in residue.atoms:

            # REMOVE TYPES IF ALREADY ASSIGNED FROM SMARTS
            for atom_type in list(PROT_ATOM_TYPES.keys()):
                atypes_array[atom.index, atypes[atom_type]] = 0

            # ADD ATOM TYPES FROM DICTIONARY
            for atom_type, atom_ids in PROT_ATOM_TYPES.items():
                atom_id = residue.resname.strip() + atom.name.strip()
                if atom_id in atom_ids:
                    atypes_array[atom.index, atypes[atom_type]] = 1

    return atypes_array


def assign_radii(mol):
    """Get the van der Waals and covalent radii for each atom in the molecule.

    Parameters
    ----------
    mol : openbabel.OBMol
        The molecule object.

    Returns
    -------
    radii_array : np.ndarray
        An array of shape (n_atoms, 2) where the first column is the van der Waals
        radius and the second column is the covalent radius.
    """
    atom_radii = np.zeros(mol.NumAtoms())
    for atom in ob.OBMolAtomIter(mol):
        atom_radii[atom.GetId()] = ob.GetVdwRad(atom.GetAtomicNum())
    return atom_radii


def find_hydrogen_bonded_atoms(mol):
    """
    Find hydrogen bonded atoms in the molecule.
    """
    # hbond_indices = np.zeros(mol.NumAtoms(), dtype=bool)
    hbond_array = np.zeros((mol.NumAtoms(), 6), dtype=int)

    for atom in ob.OBMolAtomIter(mol):
        if atom.ExplicitHydrogenCount():
            for ix, atom2 in enumerate(ob.OBAtomAtomIter(atom)):
                if atom2.GetAtomicNum() == 1:
                    # hbond_indices[atom.GetId()] = True
                    hbond_array[atom.GetId(), ix] = atom2.GetId()

    return hbond_array
