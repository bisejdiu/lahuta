"""
Placehoder for the atom types and radii.
"""

from typing import Dict

import numpy as np
from MDAnalysis.topology.tables import vdwradii as MDA_VDW_RADII  # type: ignore
from numpy.typing import NDArray
from openbabel import openbabel as ob  # type: ignore

from lahuta.config.atoms import ID_TO_TYPES, PROT_ATOM_TYPES, STANDARD_AMINO_ACIDS
from lahuta.config.smarts import AVAILABLE_ATOM_TYPES, SmartsPatternRegistry
from lahuta.types.mdanalysis import AtomGroupType
from lahuta.types.openbabel import MolType, ObSmartPatternType, OBSmartsPatternWrapper


def assign_atom_types(mol: MolType, atomgroup: AtomGroupType) -> NDArray[np.int8]:
    """
    Assign atom types to each atom in the molecule.
    Atom types are defined in `SmartsPatternRegistry`
    """
    atypes = AVAILABLE_ATOM_TYPES
    # uv = atomgroup.universe.atoms.universe

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
            for prot_atype in list(PROT_ATOM_TYPES.keys()):
                atypes_array[atom.index, atypes[prot_atype]] = 0

            # ADD ATOM TYPES FROM DICTIONARY
            for prot_atype, atom_ids in PROT_ATOM_TYPES.items():
                atom_id = residue.resname.strip() + atom.name.strip()
                if atom_id in atom_ids:
                    atypes_array[atom.index, atypes[prot_atype]] = 1

    return atypes_array


def vec_assign_atom_types(
    mol: MolType, atomgroup: AtomGroupType, ta: Dict[str, NDArray[np.str_]]
) -> NDArray[np.int8]:
    """
    Assign atom types to each atom in the molecule.
    Atom types are defined in `SmartsPatternRegistry`
    """
    # atypes = AVAILABLE_ATOM_TYPES

    # atom_id_to_type_index = {atom_id: atypes[atom_type] for atom_type, atom_ids in PROT_ATOM_TYPES_SET.items() for atom_id in atom_ids}
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

                if atom.GetResidue().GetName() not in STANDARD_AMINO_ACIDS:
                    atypes_array[atom.GetId(), atypes[atom_type.name]] = 1

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
    atom_ids: NDArray[np.int32] = np.core.defchararray.add(  # type: ignore
        np.core.defchararray.strip(resname_str),  # type: ignore
        np.core.defchararray.strip(atom_name_str),  # type: ignore
    )

    for idx, atom in enumerate(ag):
        atom_id: int = atom_ids[idx]  # type: ignore
        atom_types = ID_TO_TYPES.get(atom_id, None)  # type: ignore

        if atom_types is not None:
            for atom_type in atom_types:
                atypes_array[atom.index, atypes[atom_type]] = 1  # type: ignore

    return atypes_array


def assign_radii(mol: MolType) -> NDArray[np.float_]:
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
        atom_radii[atom.GetId()] = ob.GetVdwRad(atom.GetAtomicNum())  # type: ignore
    return atom_radii


def find_hydrogen_bonded_atoms(mda: AtomGroupType, mol: MolType) -> NDArray[np.int32]:
    """
    Find hydrogen bonded atoms in the molecule.
    """
    # hbond_indices = np.zeros(mol.NumAtoms(), dtype=bool)
    # hbond_array = np.zeros((mol.NumAtoms(), 6), dtype=int)
    # mda = luni.to("mda").atoms
    n_atoms: int = mda.universe.atoms.n_atoms
    hbond_array: NDArray[np.int32] = np.zeros((n_atoms, 6), dtype=int)
    # mol = luni.to("mol")

    # will give -1 for atoms not in the atomgroup
    max_index = np.max(mda.indices)
    atom_mapping = np.full(max_index + 1, -1)
    atom_mapping[np.arange(mda.n_atoms)] = mda.indices

    for atom in ob.OBMolAtomIter(mol):
        if atom.ExplicitHydrogenCount():
            for ix, atom2 in enumerate(ob.OBAtomAtomIter(atom)):
                if atom2.GetAtomicNum() == 1:
                    atom1_id = atom_mapping[atom.GetId()]
                    atom2_id = atom_mapping[atom2.GetId()]
                    hbond_array[atom1_id, ix] = atom2_id

    return hbond_array


def v_radii_assignment(elements: NDArray[np.str_]) -> NDArray[np.float_]:
    vdwradii: Dict[str, int] = {k.capitalize(): v for k, v in MDA_VDW_RADII.items()}

    def v_capitalize(
        array: NDArray[np.str_], mapping: Dict[str, int]
    ) -> NDArray[np.float_]:
        vfunc = np.vectorize(mapping.get)
        return vfunc(array)

    return v_capitalize(elements, vdwradii)
