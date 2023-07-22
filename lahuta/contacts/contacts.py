import numpy as np
from numpy.typing import NDArray
from openbabel import openbabel as ob  # type: ignore

from lahuta.config.atoms import METALS
from lahuta.config.defaults import CONTACTS
from lahuta.core.neighbors import NeighborPairs
from lahuta.types.openbabel import BondIterable, MolType
from lahuta.utils.array_utils import difference, np_optimized_matching_pairs

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


def get_bonded_atoms(mol: MolType) -> NDArray[np.int32]:
    """Get the bonded atoms of a molecule.

    Parameters
    ----------
    mol : openbabel.OBMol
        The molecule object.

    Returns
    -------
    bonds : np.ndarray
        An array of shape (n_bonds, 2) where each row contains the indices of the atoms
        in the bond.
    """

    def bond_iter_wrapper(mol: MolType) -> BondIterable:
        """Wrapper for the openbabel bond iterator."""
        return ob.OBMolBondIter(mol)  # type: ignore

    bonds: NDArray[np.int32] = np.zeros((mol.NumBonds(), 2), dtype=int)
    for ix, bond in enumerate(bond_iter_wrapper(mol)):
        # assert isinstance(bond, MolBond), "bond is not a MolBond"
        atom_idx1, atom_idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonds[ix, :] = (
            (atom_idx1, atom_idx2) if atom_idx1 < atom_idx2 else (atom_idx2, atom_idx1)
        )

    return bonds


def covalent_neighbors(ns: NeighborPairs) -> NeighborPairs:
    """Find neighbors that are covalently bonded.

    Parameters
    ----------
    ns : NeighborPairs
        A NeighborPairs object.

    Returns
    -------
    NeighborPairs
        A NeighborPairs object containing only covalent neighbors.
    """
    bonds = get_bonded_atoms(ns.mol)
    indices = np_optimized_matching_pairs(ns.pairs + 1, bonds)

    return ns.clone(ns.pairs[indices], ns.distances[indices])


def metalic_neighbors(
    ns: NeighborPairs, distance: float = CONTACTS["metal"]["distance"]
) -> NeighborPairs:
    """Find neighbors that form metalic contacts.

    Parameters
    ----------
    ns : NeighborPairs
        A NeighborPairs object.

    distance : float
        The distance cutoff for metalic contacts. Check the default value in
        `lahuta.config.defaults`.

    Returns
    -------
    NeighborPairs
        A NeighborPairs object containing only metalic contacts.
    """
    metal_indices = (
        ns.atoms[ns.indices].select_atoms("element " + " ".join(METALS)).indices
    )

    acceptor_metal = (
        ns.type_filter("hbond_acceptor", 1)
        .index_filter(metal_indices, 2)
        .distance_filter(distance)
    )

    metal_acceptor = (
        ns.type_filter("hbond_acceptor", 2)
        .index_filter(metal_indices, 1)
        .distance_filter(distance)
    )

    return acceptor_metal + metal_acceptor


def carbonyl_neighbors(
    ns: NeighborPairs, distance: float = CONTACTS["carbonyl"]["distance"]
) -> NeighborPairs:
    """Find neighbors that form carbonyl contacts.

    Parameters
    ----------
    ns : NeighborPairs
        A NeighborPairs object.

    distance : float
        The distance cutoff for carbonyl contacts.

    Returns
    -------
    NeighborPairs
        A NeighborPairs object containing only carbonyl contacts.
    """

    contacts_atom12 = (
        ns.type_filter("carbonyl_oxygen", 1)
        .type_filter("carbonyl_carbon", 2)
        .distance_filter(distance)
    )

    contacts_atom21 = (
        ns.type_filter("carbonyl_carbon", 1)
        .type_filter("carbonyl_oxygen", 2)
        .distance_filter(distance)
    )

    return contacts_atom12 + contacts_atom21


def ionic_neighbors(
    ns: NeighborPairs, distance: float = CONTACTS["ionic"]["distance"]
) -> NeighborPairs:
    """Find neighbors that form ionic contacts.

    Parameters
    ----------
    ns : NeighborPairs
        A NeighborPairs object.

    distance : float
        The distance cutoff for ionic contacts. Check the default value in
        `lahuta.config.defaults`.

    Returns
    -------
    NeighborPairs
        A NeighborPairs object containing only ionic contacts.

    """
    contacts_atom12 = (
        ns.type_filter("pos_ionisable", 1)
        .type_filter("neg_ionisable", 2)
        .distance_filter(distance)
    )

    contacts_atom21 = (
        ns.type_filter("neg_ionisable", 1)
        .type_filter("pos_ionisable", 2)
        .distance_filter(distance)
    )

    return contacts_atom12 + contacts_atom21


def aromatic_neighbors(
    ns: NeighborPairs, distance: float = CONTACTS["aromatic"]["distance"]
) -> NeighborPairs:
    """Find neighbors that form aromatic contacts.

    Parameters
    ----------
    ns : NeighborPairs
        A NeighborPairs object.

    distance : float
        The distance cutoff for aromatic contacts. Check the default value in
        `lahuta.config.defaults`.

    Returns
    -------
    NeighborPairs
        A NeighborPairs object containing only aromatic contacts.
    """
    return (
        ns.type_filter("aromatic", 1)
        .type_filter("aromatic", 2)
        .distance_filter(distance)
    )


def hydrophobic_neighbors(
    ns: NeighborPairs, distance: float = CONTACTS["hydrophobic"]["distance"]
) -> NeighborPairs:
    """Find neighbors that form hydrophobic contacts.

    Parameters
    ----------
    ns : NeighborPairs
        A NeighborPairs object.

    distance : float
        The distance cutoff for hydrophobic contacts. Check the default value in
        `lahuta.config.defaults`.

    Returns
    -------
    NeighborPairs
        A NeighborPairs object containing only hydrophobic contacts.
    """
    return (
        ns.type_filter("hydrophobe", 1)
        .type_filter("hydrophobe", 2)
        .distance_filter(distance)
    )


def vdw_neighbors(
    ns: NeighborPairs, vdw_comp_factor: float = 0.1, remove_clashes: bool = True
) -> NeighborPairs:
    """Find neighbors that form vdw contacts.

    Parameters
    ----------
    ns : NeighborPairs
        A NeighborPairs object.

    vdw_comp_factor : float
        The factor by which the vdw radii are increased to form the vdw contact.

    remove_clashes : bool
        If True, remove clashes from the NeighborPairs object. Clashes are defined as contact
        distances less than the sum of the vdw radii of the two atoms.

    Returns
    -------
    NeighborPairs
        A NeighborPairs object containing only vdw contacts.
    """

    vdw_radii = ns.atoms.vdw_radii[ns.pairs[:, 0]] + ns.atoms.vdw_radii[ns.pairs[:, 1]]

    distance_mask = ns.distances <= vdw_radii + vdw_comp_factor
    vdw_comp_pairs = ns.pairs[distance_mask]
    vdw_distances = ns.distances[distance_mask]

    if not remove_clashes:
        return ns.clone(vdw_comp_pairs, vdw_distances)  # TODO: check if this is correct
        # return NeighborPairs(ns.luni, vdw_comp_pairs, vdw_distances)

    vdw_clash_pairs = ns.pairs[ns.distances < vdw_radii]
    no_clash_indices = difference(vdw_comp_pairs, vdw_clash_pairs)

    return ns.clone(
        vdw_comp_pairs[no_clash_indices], vdw_distances[no_clash_indices]
    )  # TODO: check if this is correct
    # return NeighborPairs(
    #     ns.luni, vdw_comp_pairs[no_clash_indices], vdw_distances[no_clash_indices]
    # )


def hbond_neighbors(ns: NeighborPairs) -> NeighborPairs:
    """Find neighbors that form hydrogen bonds.

    Parameters
    ----------
    ns : NeighborPairs
        A NeighborPairs object.

    Returns
    -------
    NeighborPairs
        A NeighborPairs object containing only hydrogen bonds.
    """

    hbond_atom12 = (
        ns.type_filter("hbond_donor", 1)
        .type_filter("hbond_acceptor", 2)
        .hbond_distance_filter(partner=2)
        .hbond_angle_filter(partner=1)
    )

    hbond_atom21 = (
        ns.type_filter("hbond_donor", 2)
        .type_filter("hbond_acceptor", 1)
        .hbond_distance_filter(partner=1)
        .hbond_angle_filter(partner=2)
    )

    return hbond_atom12 + hbond_atom21


def weak_hbond_neighbors(ns: NeighborPairs) -> NeighborPairs:
    """Find neighbors that form weak hydrogen bonds.

    Parameters
    ----------
    ns : NeighborPairs
        A NeighborPairs object.

    Returns
    -------
    NeighborPairs
        A NeighborPairs object containing only weak hydrogen bonds.
    """

    hbond_atom12 = (
        ns.type_filter("hbond_acceptor", 1)
        .type_filter("weak_hbond_donor", 2)
        .hbond_distance_filter(partner=1)
        .hbond_angle_filter(partner=2, weak=True)
    )

    hbond_atom21 = (
        ns.type_filter("hbond_acceptor", 2)
        .type_filter("weak_hbond_donor", 1)
        .hbond_distance_filter(partner=2)
        .hbond_angle_filter(partner=1, weak=True)
    )

    return hbond_atom12 + hbond_atom21


def polar_hbond_neighbors(
    ns: NeighborPairs, distance: float = CONTACTS["hbond"]["polar distance"]
) -> NeighborPairs:
    """Find neighbors that form polar hydrogen bonds.

    Parameters
    ----------
    ns : NeighborPairs
        A NeighborPairs object.

    distance : float
        The distance cutoff for polar hydrogen bonds. Check the default value in
        `lahuta.config.defaults`.

    Returns
    -------
    NeighborPairs
        A NeighborPairs object containing only polar hydrogen bonds.
    """

    hbond_atom12 = (
        ns.type_filter("hbond_donor", 1)
        .type_filter("hbond_acceptor", 2)
        .distance_filter(distance)
    )

    hbond_atom21 = (
        ns.type_filter("hbond_donor", 2)
        .type_filter("hbond_acceptor", 1)
        .distance_filter(distance)
    )

    return hbond_atom12 + hbond_atom21


def weak_polar_hbond_neighbors(
    ns: NeighborPairs, distance: float = CONTACTS["weak hbond"]["weak polar distance"]
) -> NeighborPairs:
    """Find neighbors that form weak polar hydrogen bonds.

    Parameters
    ----------
    ns : NeighborPairs
        A NeighborPairs object.

    distance : float
        The distance cutoff for weak polar hydrogen bonds. Check the default value in
        `lahuta.config.defaults`.

    Returns
    -------
    NeighborPairs
        A NeighborPairs object containing only weak polar hydrogen bonds.
    """

    hbond_atom12 = (
        ns.type_filter("hbond_acceptor", 1)
        .type_filter("weak_hbond_donor", 2)
        .distance_filter(distance)
    )

    hbond_atom21 = (
        ns.type_filter("hbond_acceptor", 2)
        .type_filter("weak_hbond_donor", 1)
        .distance_filter(distance)
    )

    return hbond_atom12 + hbond_atom21
