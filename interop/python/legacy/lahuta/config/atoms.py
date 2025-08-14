"""Provides definitions for various sets of atoms and atom types used
in chemical and biochemical contexts. Additionally, it includes a helper function 
for remapping atom types for easier lookups.

Variables:
    ```
    HALOGENS (set): A set of standard halogens.
    MAINCHAIN_ATOMS (set): A set of atoms usually found in the main chain of a protein.
    STANDARD_NUCLEOTIDES (set): A set of standard nucleotides in RNA and DNA.
    METALS (set): A set of metal atoms, sourced from `atom_types.py`.
    PROTEIN_RESIDUES (set): A set of standard amino acids, sourced from `atom_types.py`.
    PROT_ATOM_TYPES (dict): A dictionary mapping atom types to sets of atom names (strings),
        as defined in `atom_types.py`.
    ID_TO_TYPES (dict): A dictionary mapping atom names (strings) to sets of atom types,
        as defined in `atom_types.py`.
    ```

The module is expected to be used for handling and processing biochemical data structures, 
with specific focuses on proteins and nucleic acids.
"""


import lahuta.config.atom_types as at

HALOGENS = {"F", "CL", "BR", "I", "AT"}
"""
Type: `set[str]`: A set of standard halogens.
"""

MAINCHAIN_ATOMS = {"N", "C", "CA", "O", "OXT"}
"""
Type: `set[str]`: A set of atoms usually found in the main chain of a protein.
"""


STANDARD_NUCLEOTIDES = {"A", "C", "G", "I", "U", "DA", "DC", "DG", "DI", "DT", "DU", "N"}
"""
Type: `set[str]`: A set of standard nucleotides in RNA and DNA.
"""


METALS = at.METALS
"""
Type: `set[str]`: A set of metal atoms, sourced from `lahuta.config.atom_types.py`.
"""

PROTEIN_RESIDUES = at.PROTEIN_RESIDUES
"""
Type: `set[str]`: A set of standard amino acids. See the 
[source](atom_types.md#lahuta.config.atom_types.PROTEIN_RESIDUES).
"""


PROT_ATOM_TYPES: dict[str, set[str]] = {
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
"""
Type: `dict[str, set[str]]`: A dictionary mapping atom types to sets of atom names (strings).

!!! tip "Definitions"

    - `hbond_acceptor` -> [`_HA_ATOM_TYPES`](atom_types.md#lahuta.config.atom_types._HA_ATOM_TYPES)
    - `hbond_donor` -> [`_HD_ATOM_TYPES`](atom_types.md#lahuta.config.atom_types._HD_ATOM_TYPES)
    - `xbond_acceptor` -> [`_XA_ATOM_TYPES`](atom_types.md#lahuta.config.atom_types._XA_ATOM_TYPES)
    - `xbond_donor` -> [`_XD_ATOM_TYPES`](atom_types.md#lahuta.config.atom_types._XD_ATOM_TYPES)
    - `weak_hbond_acceptor` -> [`_WHA_ATOM_TYPES`](atom_types.md#lahuta.config.atom_types._WHA_ATOM_TYPES)
    - `weak_hbond_donor` -> [`_WHD_ATOM_TYPES`](atom_types.md#lahuta.config.atom_types._WHD_ATOM_TYPES)
    - `pos_ionisable` -> \
        [`_POS_IONISABLE_ATOM_TYPES`](atom_types.md#lahuta.config.atom_types._POS_IONISABLE_ATOM_TYPES)
    - `neg_ionisable` -> \
        [`_NEG_IONISABLE_ATOM_TYPES`](atom_types.md#lahuta.config.atom_types._NEG_IONISABLE_ATOM_TYPES)
    - `hydrophobe` -> [`_HYDROPHOBE_ATOM_TYPES`](atom_types.md#lahuta.config.atom_types._HYDROPHOBE_ATOM_TYPES)
    - `carbonyl_oxygen` -> \
        [`_CARBONYL_OXYGEN_ATOM_TYPES`](atom_types.md#lahuta.config.atom_types._CARBONYL_OXYGEN_ATOM_TYPES)
    - `carbonyl_carbon` -> \
        [`_CARBONYL_CARBON_ATOM_TYPES`](atom_types.md#lahuta.config.atom_types._CARBONYL_CARBON_ATOM_TYPES)
    - `aromatic` -> [`_AROMATIC_ATOM_TYPES`](atom_types.md#lahuta.config.atom_types._AROMATIC_ATOM_TYPES)

"""


def remap_prot_atom_types(prot_atom_types: dict[str, set[str]]) -> dict[str, set[str]]:
    """Convert the PROT_ATOM_TYPES dictionary to a dictionary mapping atom ids to atom types.

    Args:
        prot_atom_types (dict): A dictionary mapping atom types to atom ids.

    Returns:
        dict: A dictionary mapping atom ids to atom types.
    """
    id_to_types: dict[str, set[str]] = {}
    for atom_type, atom_set in prot_atom_types.items():
        for atom_id in atom_set:
            if atom_id not in id_to_types:
                id_to_types[atom_id] = set()
            id_to_types[atom_id].add(atom_type)
    return id_to_types


ID_TO_TYPES = remap_prot_atom_types(PROT_ATOM_TYPES)
"""
Type: `dict[str, set[str]]`: A dictionary mapping atom names (strings) to sets of atom types.

!!! tip "Definitions"
    See [`PROT_ATOM_TYPES`](#lahuta.config.atoms.PROT_ATOM_TYPES) for the definitions of atom types.

"""
