"""Implements the core logic for calculating different types of atomic contacts in proteins. It provides
a variety of functions, each corresponding to a specific type of atomic contact including covalent, metallic, 
carbonyl, ionic, aromatic, hydrophobic, van der Waals, and different kinds of hydrogen bonds. 

Functions:
    covalent_neighbors(ns: NeighborPairs) -> NeighborPairs: Computes covalent contacts.
    metalic_neighbors(ns: NeighborPairs) -> NeighborPairs: Computes metallic contacts.
    carbonyl_neighbors(ns: NeighborPairs) -> NeighborPairs: Computes carbonyl contacts.
    ionic_neighbors(ns: NeighborPairs) -> NeighborPairs: Computes ionic contacts.
    aromatic_neighbors(ns: NeighborPairs) -> NeighborPairs: Computes aromatic contacts.
    hydrophobic_neighbors(ns: NeighborPairs) -> NeighborPairs: Computes hydrophobic contacts.
    vdw_neighbors(ns: NeighborPairs) -> NeighborPairs: Computes van der Waals contacts.
    hbond_neighbors(ns: NeighborPairs) -> NeighborPairs: Computes hydrogen bonds.
    weak_hbond_neighbors(ns: NeighborPairs) -> NeighborPairs: Computes weak hydrogen bonds.
    polar_hbond_neighbors(ns: NeighborPairs) -> NeighborPairs: Computes polar hydrogen bonds.
    weak_polar_hbond_neighbors(ns: NeighborPairs) -> NeighborPairs: Computes weak polar hydrogen bonds.

Each function takes a `NeighborPairs` object, which represents precomputed neighbor relationships between atoms, 
and an optional distance cutoff parameter for certain types of contacts. Each function returns a new `NeighborPairs` 
object that represents the specific contacts computed by that function.

Notes:
    Each atomic contact type is computed based on a set of conditions, primarily involving the types of the two 
    interacting atoms and the distance between them. For example, carbonyl contacts are calculated by filtering 
    for neighbor pairs where one atom is a carbonyl oxygen and the other is a carbonyl carbon, and the distance 
    between them is less than a predefined cutoff.
    
    Most contact calculation functions also accept a `distance` argument specifying the cutoff distance for 
    considering two atoms as neighbors in the contact type calculation. These cutoffs are defined in the 
    `lahuta.config.defaults` module and can be adjusted as per user requirements.

"""

from typing import Optional

from lahuta.config.atoms import METALS
from lahuta.config.defaults import CONTACTS
from lahuta.core._hbond_handler import HBondHandler
from lahuta.core.neighbors import NeighborPairs
from lahuta.utils.array_utils import difference, find_shared_pairs
from lahuta.utils.ob import get_bonded_atoms

__all__ = [
    "covalent_neighbors",
    "metalic_neighbors",
    "carbonyl_neighbors",
    "ionic_neighbors",
    "aromatic_neighbors",
    "hydrophobic_neighbors",
    "vdw_neighbors",
    "hbond_neighbors",
    "weak_hbond_neighbors",
    "polar_hbond_neighbors",
    "weak_polar_hbond_neighbors",
]


def covalent_neighbors(ns: NeighborPairs) -> NeighborPairs:
    """Handle the computation of covalent contacts in a molecular system.

    Covalent contacts are interactions based on covalent bonds between atoms in a molecular system.
    This class, a derivative of the `ContactAnalysis` base class, overrides the `compute` method
    to provide functionality specifically for covalent contact computation.

    Covalent contacts refer to the interactions between atoms that share an electron pair, forming a covalent bond.
    We use the OpenBabel library to identify covalent bonds in the structure.

    !!! tip "Definition"
        Two atoms are considered to form a covalent contact if they are covalently bonded according to the molecular
        structure information obtained from OpenBabel.

    Args:
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.

    Returns:
        (NeighborPairs): A NeighborPairs object containing only covalent contacts.
    """
    bonds = get_bonded_atoms(ns.luni.to("mol"))
    indices = find_shared_pairs(ns.pairs + 1, bonds)

    return ns.new(ns.pairs[indices], ns.distances[indices])


def metalic_neighbors(ns: NeighborPairs, distance: float = CONTACTS["metal"]["distance"]) -> NeighborPairs:
    """Handle the computation of metallic contacts in a molecular system.

    Metallic contacts are interactions involving metal atoms in a molecular system.
    This class, a derivative of the `ContactAnalysis` base class, overrides the `compute`
    method to provide functionality specifically for metallic contact computation.

    Metallic contacts are interactions between metal ions and atoms that can act as hydrogen bond acceptors.
    These contacts play significant roles in the structure and function of many proteins, especially metalloproteins.
    Metal ions can form coordination bonds with electron-rich atoms (like oxygen, nitrogen, or sulfur),
    contributing to the structural stability and sometimes the catalytic activity of these proteins.

    !!! tip "Definition"
        1. The contact involves a metal ion and an atom that is a hydrogen bond acceptor.
        2. The distance between the metal ion and the hydrogen bond acceptor \
            does not exceed the predefined distance cutoff.

    Args:
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.
        distance (float): The maximum distance to consider for a metallic contact. This value is retrieved
            from the 'metal' entry of the global CONTACTS dictionary.

    Returns:
        (NeighborPairs): A NeighborPairs object containing only metallic contacts.
    """
    metal_indices = ns.atoms.select_atoms("element " + " ".join(METALS)).indices

    acceptor_metal = ns.type_filter("hbond_acceptor", 1).index_filter(metal_indices, 2).distance_filter(distance)

    metal_acceptor = ns.type_filter("hbond_acceptor", 2).index_filter(metal_indices, 1).distance_filter(distance)

    return acceptor_metal + metal_acceptor


def carbonyl_neighbors(ns: NeighborPairs, distance: float = CONTACTS["carbonyl"]["distance"]) -> NeighborPairs:
    """Handle the computation of carbonyl contacts in a molecular system.

    Carbonyl contacts involve the interaction between a carbonyl oxygen atom (O) and a carbonyl carbon atom (C)
    from a carbonyl functional group (C=O) in the context of protein-ligand structures or protein-protein structures.

    In a carbonyl group, the carbon atom has a double bond with the oxygen atom. This arrangement results in a polar
    bond with the oxygen atom carrying a partial negative charge and the carbon atom a partial positive charge.
    This polarity can lead to interactions with other polar or charged atoms.

    !!! tip "Definition"
        1. One atom is a carbonyl oxygen atom (O=C).
        2. The second atom is a carbonyl carbon atom (O=C).
        3. The distance between these two atoms does not exceed a defined distance cutoff.

    The directionality of the contact is not considered, meaning that an Oxygen to Carbon contact is equivalent
    to a Carbon to Oxygen contact.

    Args:
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.
        distance (float): The maximum distance to consider for a carbonyl contact. This value is retrieved

    Returns:
        (NeighborPairs): A NeighborPairs object containing only carbonyl contacts.
    """
    contacts_atom12 = ns.type_filter("carbonyl_oxygen", 1).type_filter("carbonyl_carbon", 2).distance_filter(distance)

    contacts_atom21 = ns.type_filter("carbonyl_carbon", 1).type_filter("carbonyl_oxygen", 2).distance_filter(distance)

    return contacts_atom12 + contacts_atom21


def ionic_neighbors(ns: NeighborPairs, distance: float = CONTACTS["ionic"]["distance"]) -> NeighborPairs:
    """Handle the computation of ionic contacts in a molecular system.

    Ionic contacts refer to the interactions between positively and negatively ionizable atoms,
    forming one of the primary types of electrostatic interactions.

    !!! tip "Definition"
        1. One atom must be positively ionizable.
        2. The other atom must be negatively ionizable.
        3. The distance between these two atoms does not exceed a defined distance cutoff.

    These criteria apply regardless of the order of the atoms in the pair, meaning a positively ionizable
    to negatively ionizable contact is considered equivalent to a negatively ionizable to positively ionizable contact.

    Args:
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.
        distance (float): The maximum distance to consider for an ionic contact. This value is retrieved from

    Returns:
        (NeighborPairs): A NeighborPairs object containing only ionic contacts.

    """
    contacts_atom12 = ns.type_filter("pos_ionisable", 1).type_filter("neg_ionisable", 2).distance_filter(distance)

    contacts_atom21 = ns.type_filter("neg_ionisable", 1).type_filter("pos_ionisable", 2).distance_filter(distance)

    return contacts_atom12 + contacts_atom21


def aromatic_neighbors(ns: NeighborPairs, distance: float = CONTACTS["aromatic"]["distance"]) -> NeighborPairs:
    """Handle the computation of aromatic contacts in a molecular structure.

    Aromatic contacts are computed based on the interactions between atoms in
    aromatic rings found in proteins and ligands. Aromatic interactions,
    commonly found in biological systems, are characterized by π-stacking
    (face-to-face), T-shaped or edge-to-face configurations, and cation-π
    interactions.

    !!! tip "Definition"
        1. Both atoms belong to an aromatic ring.
        2. The distance between these two atoms does not exceed a defined distance cutoff.

    Args:
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.
        distance (float): The maximum distance to consider for contact. Check the default value in
            `lahuta.config.defaults`.

    Returns:
        (NeighborPairs): A NeighborPairs object containing only aromatic contacts.
    """
    return ns.type_filter("aromatic", 1).type_filter("aromatic", 2).distance_filter(distance)


def hydrophobic_neighbors(ns: NeighborPairs, distance: float = CONTACTS["hydrophobic"]["distance"]) -> NeighborPairs:
    """Handle the computation of hydrophobic contacts in a molecular system.

    Hydrophobic contacts are interactions between hydrophobic atoms in a molecular system.
    This class, a derivative of the `ContactAnalysis` base class, overrides the `compute`
    method to provide functionality specifically for hydrophobic contact computation.

    Hydrophobic contacts are interactions between hydrophobic (non-polar) atoms.
    Hydrophobic interactions occur due to the tendency of hydrophobic molecules to aggregate together
    in an aqueous environment, minimizing their exposure to water molecules.

    !!! tip "Definition"
        1. Both atoms must be hydrophobic.
        2. The distance between these two atoms does not exceed a defined distance cutoff.

    Args:
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.
        distance (float): The maximum distance to consider for a hydrophobic contact.

    Returns:
        (NeighborPairs): A NeighborPairs object containing only hydrophobic contacts.
    """
    return ns.type_filter("hydrophobe", 1).type_filter("hydrophobe", 2).distance_filter(distance)


def vdw_neighbors(ns: NeighborPairs, vdw_comp_factor: float = 0.1, remove_clashes: bool = True) -> NeighborPairs:
    """Handle the computation of Van der Waals (VdW) contacts in a molecular system.

    Van der Waals (VdW) contacts are determined based on the interactions between atoms that come
    within their combined van der Waals radii, increased by a compensation factor.
    Van der Waals interactions occur due to induced polarization of atoms and can play a critical role
    in the stability of biological structures and in molecular recognition processes.

    !!! tip "Definition"
        1. The distance between two atoms does not exceed the sum of their van der Waals radii,
        increased by a defined compensation factor.
        2. Optionally, if the distance between two atoms is less than the sum of their van der Waals radii
        (defined as a clash), the contact can be excluded.

    Args:
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.
        vdw_comp_factor (float): The factor by which Van der Waals radii are multiplied during contact computation.
            This factor is used in the F.vdw_neighbors function to scale the VdW radii. Default value is 0.1.
        remove_clashes (bool): Flag indicating whether clashes (interactions with distances smaller than
            the combined Van der Waals radii of two atoms) should be removed. Default value is True.

    Returns:
        (NeighborPairs): A NeighborPairs object containing only Van der Waals contacts.
    """
    vdw_radii = ns.atoms.vdw_radii[ns.pairs[:, 0]] + ns.atoms.vdw_radii[ns.pairs[:, 1]]

    distance_mask = ns.distances <= vdw_radii + vdw_comp_factor
    vdw_comp_pairs = ns.pairs[distance_mask]
    vdw_distances = ns.distances[distance_mask]

    if not remove_clashes:
        return ns.new(vdw_comp_pairs, vdw_distances)  # TODO @bisejdiu: check if this is correct

    vdw_clash_pairs = ns.pairs[ns.distances < vdw_radii]
    no_clash_indices = difference(vdw_comp_pairs, vdw_clash_pairs)

    return ns.new(vdw_comp_pairs[no_clash_indices], vdw_distances[no_clash_indices])


def hbond_neighbors(ns: NeighborPairs, _handler: Optional[HBondHandler] = None) -> NeighborPairs:
    """Handle the computation of hydrogen bond (hbond) contacts in a molecular system.

    Hydrogen bonds are pivotal non-covalent interactions that significantly influence the structure, stability,
    and dynamics of biomolecules such as proteins and nucleic acids. A hydrogen bond forms when a hydrogen atom
    covalently bonded to a highly electronegative atom (such as oxygen or nitrogen) interacts with another
    electronegative atom from a different group.

    !!! tip "Definition"
        1. It involves a hydrogen bond donor atom and a hydrogen bond acceptor atom.
        2. The distance between the hydrogen atom and the acceptor atom does not exceed a predefined cutoff distance.
        3. The angle formed by the donor, hydrogen, and acceptor atoms falls within a predefined range.
        This accounts for the directional nature of hydrogen bonds.

    Args:
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.

    Returns:
        (NeighborPairs): A NeighborPairs object containing only hbond contacts.
    """
    handler = _handler or HBondHandler(ns)

    hbond_atom12 = (
        ns.type_filter("hbond_donor", 1)
        .type_filter("hbond_acceptor", 2)
    )
    hbond_atom12 = handler.hbond_distance_filter(hbond_atom12, partner=2)
    hbond_atom12 = handler.hbond_angle_filter(hbond_atom12, partner=1)

    hbond_atom21 = (
        ns.type_filter("hbond_donor", 2)
        .type_filter("hbond_acceptor", 1)
    )
    hbond_atom21 = handler.hbond_distance_filter(hbond_atom21, partner=1)
    hbond_atom21 = handler.hbond_angle_filter(hbond_atom21, partner=2)

    return hbond_atom12 + hbond_atom21


def weak_hbond_neighbors(ns: NeighborPairs, _handler: Optional[HBondHandler] = None) -> NeighborPairs:
    """Handle the computation of weak hydrogen bond (weak hbond) contacts in a molecular system.

    Weak hydrogen bonds are a type of non-covalent interactions that, despite their reduced strength
    compared to regular hydrogen bonds, still play important roles in biomolecular structures and functions.
    These bonds form when a hydrogen atom, covalently linked to a weakly electronegative atom (such as carbon),
    interacts with an electronegative acceptor atom from a different group. The difference in electronegativity
    between the donor and acceptor atoms is what makes these bonds relatively weaker than conventional hydrogen bonds.

    !!! tip "Definition"
        1. A weak hydrogen bond involves a weak hydrogen bond donor atom and a hydrogen bond acceptor atom.
        2. The distance between the hydrogen atom and the acceptor atom does not surpass a certain cutoff distance.
        3. The angle formed by the donor, hydrogen, and acceptor atoms is within a predefined range.
        This accounts for the directional nature of hydrogen bonds.

    Args:
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.

    Returns:
        (NeighborPairs): A NeighborPairs object containing only weak hydrogen bonds.
    """
    handler = _handler or HBondHandler(ns)
    hbond_atom12 = (
        ns.type_filter("hbond_acceptor", 1)
        .type_filter("weak_hbond_donor", 2)
    )
    hbond_atom12 = handler.hbond_distance_filter(hbond_atom12, partner=1)
    hbond_atom12 = handler.hbond_angle_filter(hbond_atom12, partner=2, weak=True)

    hbond_atom21 = (
        ns.type_filter("hbond_acceptor", 2)
        .type_filter("weak_hbond_donor", 1)
    )
    hbond_atom21 = handler.hbond_distance_filter(hbond_atom21, partner=2)
    hbond_atom21 = handler.hbond_angle_filter(hbond_atom21, partner=1, weak=True)

    return hbond_atom12 + hbond_atom21


def polar_hbond_neighbors(ns: NeighborPairs, distance: float = CONTACTS["hbond"]["polar distance"]) -> NeighborPairs:
    """Handle the computation of polar hydrogen bond (polar hbond) contacts in a molecular system.

    Polar hydrogen bonds involve a hydrogen atom covalently bonded to a polar atom (hydrogen bond donor),
    forming an interaction with another polar atom from a different group (hydrogen bond acceptor). In contrast to
    conventional hydrogen bonds, for polar hydrogen bonds, we do not consider the angle formed by the donor, hydrogen,
    and acceptor atoms.

    !!! tip "Definition"
        1. The presence of a hydrogen bond donor and a hydrogen bond acceptor.
        2. The distance between the hydrogen bond donor and the hydrogen bond acceptor does
        not exceed a specific cutoff distance.

    Args:
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.
        distance (float): The maximum distance to consider for a polar hydrogen bond. Check the default value in
            `lahuta.config.defaults`.

    Returns:
        (NeighborPairs): A NeighborPairs object containing only polar hydrogen bonds.
    """
    hbond_atom12 = ns.type_filter("hbond_donor", 1).type_filter("hbond_acceptor", 2).distance_filter(distance)

    hbond_atom21 = ns.type_filter("hbond_donor", 2).type_filter("hbond_acceptor", 1).distance_filter(distance)

    return hbond_atom12 + hbond_atom21


def weak_polar_hbond_neighbors(
    ns: NeighborPairs, distance: float = CONTACTS["weak hbond"]["weak polar distance"]
) -> NeighborPairs:
    """Handle the computation of weak polar hydrogen bond (weak polar hbond) contacts in a molecular system.

    Weak polar hydrogen bonds rely on a weak hydrogen bond donor, and we also do not consider the angle
    formed by the donor, hydrogen, and acceptor atoms.

    !!! tip "Definition"
        1. The presence of a hydrogen bond acceptor and a weak hydrogen bond donor.
        2. The distance between the hydrogen bond acceptor and the weak hydrogen bond donor
        should not exceed a certain cutoff distance.

    Args:
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.
        distance (float): The maximum distance to consider for a weak polar hydrogen bond. Check the default value in
            `lahuta.config.defaults`.

    Returns:
        (NeighborPairs): A NeighborPairs object containing only weak polar hydrogen bonds.ff
    """
    hbond_atom12 = ns.type_filter("hbond_acceptor", 1).type_filter("weak_hbond_donor", 2).distance_filter(distance)

    hbond_atom21 = ns.type_filter("hbond_acceptor", 2).type_filter("weak_hbond_donor", 1).distance_filter(distance)

    return hbond_atom12 + hbond_atom21
