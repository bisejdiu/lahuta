from typing import Dict, Set

import lahuta.config._atom_type_strings as at

HALOGENS = set(["F", "CL", "BR", "I", "AT"])

MAINCHAIN_ATOMS = set(["N", "C", "CA", "O", "OXT"])

STANDARD_NUCLEOTIDES = set(
    ["A", "C", "G", "I", "U", "DA", "DC", "DG", "DI", "DT", "DU", "N"]
)

METALS = at.METALS
STANDARD_AMINO_ACIDS = at.STANDARD_AMINO_ACIDS


PROT_ATOM_TYPES: Dict[str, Set[str]] = {
    "hbond_acceptor": at.HBOND_ACCEPTORS,
    "hbond_donor": at.HBOND_DONORS,
    "xbond_acceptor": at.XBOND_ACCEPTORS,
    "xbond_donor": set(),
    "weak_hbond_acceptor": at.WEAK_HBOND_ACCEPTORS,
    "weak_hbond_donor": at.WEAK_HBOND_DONORS,
    "pos_ionisable": at.POS_IONISABLE,
    "neg_ionisable": at.NEG_IONISABLE,
    "hydrophobe": at.HYDROPHOBES,
    "carbonyl_oxygen": at.CARBONYL_OXYGENS,
    "carbonyl_carbon": at.CARBONYL_CARBONS,
    "aromatic": at.AROMATIC,
}


def remap_prot_atom_types(prot_atom_types):
    """Convert the PROT_ATOM_TYPES dictionary to a dictionary mapping atom ids to atom types.

    Args:
        prot_atom_types (dict): A dictionary mapping atom types to atom ids.

    Returns:
        dict: A dictionary mapping atom ids to atom types.
    """
    id_to_types = {}
    for atom_type, atom_set in prot_atom_types.items():
        for atom_id in atom_set:
            if atom_id not in id_to_types:
                id_to_types[atom_id] = set()
            id_to_types[atom_id].add(atom_type)
    return id_to_types


ID_TO_TYPES = remap_prot_atom_types(PROT_ATOM_TYPES)
