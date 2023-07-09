"""
Placehoder for the atom types and radii.
"""
import numpy as np
from openbabel import openbabel as ob

from lahuta.config.atoms import ID_TO_TYPES

from ..config.atoms import PROT_ATOM_TYPES
from ..config.residues import STANDARD_AMINO_ACIDS
from ..config.smarts import ATOM_TYPES


def assign_atom_types(mol, atomgroup):
    """
    Assign atom types to each atom in the molecule.
    Atom types are defined in `ATOM_TYPES`
    """
    atypes = {
        "hbond_acceptor": 0,
        "pos_ionisable": 1,
        "carbonyl_oxygen": 2,
        "weak_hbond_donor": 3,
        "carbonyl_carbon": 4,
        "weak hbond_acceptor": 5,
        "hbond_donor": 6,
        "neg_ionisable": 7,
        "aromatic": 8,
        "xbond_acceptor": 9,
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
        atypes_array[atom.index, atypes["hbond_acceptor"]] = 1
        atypes_array[atom.index, atypes["hbond_donor"]] = 1

    # OVERRIDE PROTEIN ATOM TYPING FROM DICTIONARY
    for residue in atomgroup.select_atoms(
        "resname " + " ".join(STANDARD_AMINO_ACIDS)
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


def vec_assign_atom_types(mol, atomgroup, ta):
    """
    Assign atom types to each atom in the molecule.
    Atom types are defined in `ATOM_TYPES`
    """
    atypes = {
        "hbond_acceptor": 0,
        "pos_ionisable": 1,
        "carbonyl_oxygen": 2,
        "weak_hbond_donor": 3,
        "carbonyl_carbon": 4,
        "weak hbond_acceptor": 5,
        "hbond_donor": 6,
        "neg_ionisable": 7,
        "aromatic": 8,
        "xbond_acceptor": 9,
        "hydrophobe": 10,
    }

    # atom_id_to_type_index = {atom_id: atypes[atom_type] for atom_type, atom_ids in PROT_ATOM_TYPES_SET.items() for atom_id in atom_ids}
    atypes = {x: i for i, x in enumerate(list(PROT_ATOM_TYPES.keys()))}

    atypes_array = np.zeros((mol.NumAtoms(), len(atypes)))
    for atom_type, smartsdict in ATOM_TYPES.items():
        for smarts in smartsdict.values():
            ob_smart = ob.OBSmartsPattern()
            ob_smart.Init(str(smarts))
            ob_smart.Match(mol)

            matches = [x[0] for x in ob_smart.GetMapList()]
            for match in matches:
                atom = mol.GetAtom(match)

                if atom.GetResidue().GetName() not in STANDARD_AMINO_ACIDS:
                    atypes_array[atom.GetId(), atypes[atom_type]] = 1

    # ALL WATER MOLECULES ARE HYDROGEN BOND DONORS AND ACCEPTORS
    for atom in atomgroup.select_atoms(
        "resname SOL HOH TIP3 TIP4 WAT W and not name H*"
    ):
        atypes_array[atom.index, atypes["hbond_acceptor"]] = 1
        atypes_array[atom.index, atypes["hbond_donor"]] = 1

    # OVERRIDE PROTEIN ATOM TYPING FROM DICTIONARY
    resname, atom_name = ta["resname"], ta["name"]

    # Convert atoms to NumPy arrays for efficient indexing
    ag = atomgroup.select_atoms("resname " + " ".join(STANDARD_AMINO_ACIDS)).atoms
    resindices = np.array([atom.resindex for atom in atomgroup])
    indices = np.array([atom.index for atom in atomgroup])

    # Convert arrays to string type
    resname_str = resname[resindices].astype(str)
    atom_name_str = atom_name[indices].astype(str)

    # Generate atom_id array by concatenating resname and atom_name arrays
    atom_ids = np.core.defchararray.add(  # type: ignore
        np.core.defchararray.strip(resname_str),  # type: ignore
        np.core.defchararray.strip(atom_name_str),  # type: ignore
    )

    for idx, atom in enumerate(ag):
        atom_id = atom_ids[idx]
        atom_types = ID_TO_TYPES.get(atom_id, None)

        if atom_types is not None:
            for atom_type in atom_types:
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
    atom_radii = np.zeros(mol.NumAtoms(), dtype=float)
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
