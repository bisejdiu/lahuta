"""Provides utility functions for assigning atom types to each atom in a molecule.
Atom types are defined in the `SmartsPatternRegistry`. 

Warning:
    As of the current version of the library, these functions are not used in the main code. 
    They are kept for comparison and testing with earlier versions of the library. 

"""

import numpy as np
from numpy.typing import NDArray
from openbabel import openbabel as ob

from lahuta._types.mdanalysis import AtomGroupType
from lahuta._types.openbabel import MolType, ObSmartPatternType, OBSmartsPatternWrapper
from lahuta.config.atoms import ID_TO_TYPES, PROT_ATOM_TYPES, PROTEIN_RESIDUES
from lahuta.config.smarts import AVAILABLE_ATOM_TYPES, SmartsPatternRegistry


def assign_atom_types(mol: MolType, atomgroup: AtomGroupType) -> NDArray[np.int8]:
    """Assign atom types to each atom in the molecule. Atom types are defined in `SmartsPatternRegistry`.

    Args:
        mol (MolType): The molecule for which to assign atom types.
        atomgroup (AtomGroupType): The atomgroup for which to assign atom types.

    Returns:
        NDArray[np.int8]: An array of shape (n_atoms, n_atom_types) where each element is either 0 or 1.
    """
    atypes = AVAILABLE_ATOM_TYPES

    atypes_array = np.zeros((mol.NumAtoms(), len(atypes)), dtype=np.int8)
    for atom_type in SmartsPatternRegistry:
        smartsdict = SmartsPatternRegistry[atom_type.name].value
        for smarts in smartsdict.values():
            ob_smart: ObSmartPatternType = OBSmartsPatternWrapper(ob.OBSmartsPattern())
            ob_smart.Init(str(smarts))
            ob_smart.Match(mol)

            matches = [x[0] for x in ob_smart.GetMapList()]
            for match in matches:
                atom = mol.GetAtom(match)

                atypes_array[atom.GetId(), atypes[atom_type.name]] = 1

    # ALL WATER MOLECULES ARE HYDROGEN BOND DONORS AND ACCEPTORS
    for atom in atomgroup.select_atoms("resname SOL HOH TIP3 TIP4 WAT W and not name H*"):
        atypes_array[atom.index, atypes["hbond_acceptor"]] = 1
        atypes_array[atom.index, atypes["hbond_donor"]] = 1

    # OVERRIDE PROTEIN ATOM TYPING FROM DICTIONARY
    for residue in atomgroup.select_atoms("resname " + " ".join(PROTEIN_RESIDUES)).residues:
        for atom in residue.atoms:
            # REMOVE TYPES IF ALREADY ASSIGNED FROM SMARTS
            for prot_atype in list(PROT_ATOM_TYPES.keys()):
                atypes_array[atom.index, atypes[prot_atype]] = 0

            # ADD ATOM TYPES FROM DICTIONARY
            for prot_atype, atom_ids in PROT_ATOM_TYPES.items():
                atom_id = residue.resname.strip() + atom.name.strip()
                if atom_id in atom_ids:
                    atypes_array[atom.index, atypes[prot_atype]] = 1

    return atypes_array


def vec_assign_atom_types(
    mol: MolType,
    atomgroup: AtomGroupType,
    ta: dict[str, NDArray[np.str_]],
) -> NDArray[np.int8]:
    """Assign atom types to each atom in the molecule. Atom types are defined in `SmartsPatternRegistry`.

    Args:
        mol (MolType): The molecule for which to assign atom types.
        atomgroup (AtomGroupType): The atomgroup for which to assign atom types.
        ta (dict[str, NDArray[np.str_]]): A dictionary containing the atom names and residue names.

    Returns:
        NDArray[np.int8]: An array of shape (n_atoms, n_atom_types) where each element is either 0 or 1.

    """
    atypes = {x: i for i, x in enumerate(list(PROT_ATOM_TYPES.keys()))}

    atypes_array = np.zeros((mol.NumAtoms(), len(atypes)), dtype=np.int8)
    for atom_type in SmartsPatternRegistry:
        smartsdict = SmartsPatternRegistry[atom_type.name].value
        for smarts in smartsdict.values():
            ob_smart: ObSmartPatternType = OBSmartsPatternWrapper(ob.OBSmartsPattern())
            ob_smart.Init(str(smarts))
            ob_smart.Match(mol)

            matches = [x[0] for x in ob_smart.GetMapList()]
            for match in matches:
                atom = mol.GetAtom(match)

                if atom.GetResidue().GetName() not in PROTEIN_RESIDUES:
                    atypes_array[atom.GetId(), atypes[atom_type.name]] = 1

    # ALL WATER MOLECULES ARE HYDROGEN BOND DONORS AND ACCEPTORS
    for atom in atomgroup.select_atoms("resname SOL HOH TIP3 TIP4 WAT W and not name H*"):
        atypes_array[atom.index, atypes["hbond_acceptor"]] = 1
        atypes_array[atom.index, atypes["hbond_donor"]] = 1

    # OVERRIDE PROTEIN ATOM TYPING FROM DICTIONARY
    resname, atom_name = ta["resname"], ta["name"]

    # Convert atoms to NumPy arrays for efficient indexing
    ag = atomgroup.select_atoms("resname " + " ".join(PROTEIN_RESIDUES)).atoms
    resindices = np.array([atom.resindex for atom in atomgroup])
    indices = np.array([atom.index for atom in atomgroup])

    # Convert arrays to string type
    resname_str = resname[resindices].astype(str)
    atom_name_str = atom_name[indices].astype(str)

    # Generate atom_id array by concatenating resname and atom_name arrays
    atom_ids: NDArray[np.str_] = np.core.defchararray.add(
        np.core.defchararray.strip(resname_str),
        np.core.defchararray.strip(atom_name_str),
    )

    for idx, atom in enumerate(ag):
        atom_id = atom_ids[idx]
        atom_types = ID_TO_TYPES.get(atom_id, None)

        if atom_types is None:
            continue

        for atom_type_x in atom_types:
            atypes_array[atom.index, atypes[atom_type_x]] = 1

    return atypes_array
