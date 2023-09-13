"""Defines atom-plane contacts, where interactions between atoms and a plane
(e.g., aromatic residues or ring structures in ligands) are computed. 

Classes:
    AtomPlaneContacts: Main class for computing atom-plane contacts.

Examples:
    ``` py
    luni = Luni(...)
    ns = luni.compute_neighbors()

    apc = AtomPlaneContacts(ns)
    dop = apc.donor_pi()    # for donor-pi contacts
    sup = apc.sulphur_pi()  # for sulphur-pi contacts
    cbp = apc.carbon_pi()   # for carbon-pi contacts
    cap = apc.cation_pi()   # for cation-pi contacts
    ```

Note:
    Functionals include a caching parameter intended to avoid re-computation 
    of shared details. However, the impact on performance isn't clear-cut; 
    thus, users are advised to test this if speed is a priority.

Warning:
    While functionals maintain consistency in the contact computation API, 
    the class-based approach is more efficient as it avoids redundant computations.
    If speed is a priority, the class-based approach is highly recommended.

"""

from typing import Callable, Optional, TypeVar

import numpy as np
from joblib import Memory
from numpy.typing import NDArray

from lahuta.config.defaults import CONTACTS
from lahuta.core.neighbors import NeighborPairs
from lahuta.lahuta_types.mdanalysis import AtomGroupType
from lahuta.utils.ob import enumerate_rings

from ._cache_funcs import (
    calc_ringnormal_pos_angle,
    calc_ringnormal_pos_angle_cached,
    compute_neighbors,
    compute_neighbors_cached,
)

memory = Memory("cachedir", verbose=0)
T = TypeVar("T")

__all__ = [
    "cation_pi",
    "carbon_pi",
    "donor_pi",
    "sulphur_pi",
    "AtomPlaneContacts",
]


DEFAULT_CONTACT_DISTS: dict[str, float] = {
    "cation_pi": CONTACTS["aromatic"]["atom_aromatic_distance"],
    "donor_pi": CONTACTS["aromatic"]["atom_aromatic_distance"],
    "sulphur_pi": CONTACTS["aromatic"]["met_sulphur_aromatic_distance"],
    "carbon_pi": CONTACTS["aromatic"]["atom_aromatic_distance"],
}


class _AtomPlaneContacts:
    """Class for computing atom-plane contacts."""

    @staticmethod
    def donor_pi(ns: NeighborPairs, angles: NDArray[np.float32], angle_cutoff: float = 30.0) -> NeighborPairs:
        """Compute the contacts between aromatic rings and the donor pi system."""
        distance = DEFAULT_CONTACT_DISTS["donor_pi"]
        return ns.numeric_filter(angles, angle_cutoff).distance_filter(distance).type_filter("hbond_donor", partner=2)

    @staticmethod
    def sulphur_pi(ns: NeighborPairs, *_: T) -> NeighborPairs:
        """Compute the contacts between aromatic rings and the sulphur pi system."""
        distance = DEFAULT_CONTACT_DISTS["sulphur_pi"]
        indices = ns.partner2.select_atoms("resname MET and element S").indices
        return ns.index_filter(indices, partner=2).distance_filter(distance)

    @staticmethod
    def carbon_pi(ns: NeighborPairs, angles: NDArray[np.float32], angle_cutoff: float = 30.0) -> NeighborPairs:
        """Compute the contacts between aromatic rings and the carbon pi system."""
        distance = DEFAULT_CONTACT_DISTS["carbon_pi"]
        return (
            ns.numeric_filter(angles, angle_cutoff)
            .index_filter(ns.partner2.select_atoms("element C").indices, partner=2)
            .distance_filter(distance)
            .type_filter("weak_hbond_donor", partner=2)
        )

    @staticmethod
    def cation_pi(ns: NeighborPairs, angles: NDArray[np.float32], angle_cutoff: float = 30.0) -> NeighborPairs:
        """Compute the contacts between aromatic rings and the cation pi system."""
        distance = DEFAULT_CONTACT_DISTS["cation_pi"]
        return ns.numeric_filter(angles, angle_cutoff).distance_filter(distance).type_filter("pos_ionisable", partner=2)


def subtract_aromatic_neighbors(
    ns: NeighborPairs, pairs: NDArray[np.int32], distances: NDArray[np.float32]
) -> NeighborPairs:
    """Subtract aromatic neighbors from the neighbor pairs."""
    cloned_neighbors = ns.clone(pairs, distances)
    return cloned_neighbors - cloned_neighbors.type_filter("aromatic", partner=2)


def compute_contacts(
    contact_fn: Callable[[NeighborPairs, NDArray[np.float32], Optional[float]], NeighborPairs],
    angle_cutoff: Optional[float],
    use_cache: bool,
) -> Callable[[NeighborPairs], NeighborPairs]:
    """Compute the contacts between aromatic rings and the specified atom plane system."""

    def wrapped(ns: NeighborPairs) -> NeighborPairs:
        """Compute the contacts between aromatic rings and the specified atom plane system."""
        mol = ns.mol
        mda: AtomGroupType = ns.luni.to("mda")
        rings = enumerate_rings(mol)

        neighbors_fn = compute_neighbors_cached if not use_cache else compute_neighbors_cached.call  # type: ignore
        result = neighbors_fn(mda.atoms.positions, rings.centers)
        pairs, distances = result[0] if use_cache else result

        neighbors = subtract_aromatic_neighbors(ns, pairs, distances)

        angles_fn = (
            calc_ringnormal_pos_angle_cached if not use_cache else calc_ringnormal_pos_angle_cached.call  # type: ignore
        )
        result = angles_fn(neighbors, mda.universe.atoms, rings.centers, rings.normals)
        angles = result[0] if use_cache else result

        return contact_fn(neighbors, angles, angle_cutoff)

    return wrapped


def create_contact_function(
    contact_type: str, angle_cutoff: Optional[float], use_cache: bool = True
) -> Callable[[NeighborPairs], NeighborPairs]:
    """Create a contact function based on the contact type."""
    func_name = f"{contact_type}"
    contact_fn = getattr(_AtomPlaneContacts, func_name)
    return compute_contacts(contact_fn, angle_cutoff, use_cache)


# user-facing functions
def cation_pi(ns: NeighborPairs, angle_cutoff: float = 30.0, cache: bool = False) -> NeighborPairs:
    """Handle the computation of cation pi contacts in a molecular system.

    Cation pi contacts are interactions between cations and the π system of atoms in a molecule.
    This class, a derivative of the `ContactAnalysis` base class, overrides the `compute`
    method to provide functionality specifically for cation pi contact computation.

    Cation-pi contacts involve interactions between a cation (a positively charged ion)
    and an electron-rich aromatic ring system. The aromatic system is capable of stabilizing the cation
    through its delocalized pi electrons.

    !!! tip "Definition"
        1. The interaction is between an aromatic ring and a positively ionizable atom.
        2. The angle between the aromatic ring plane and the vector connecting the center of the aromatic ring
        and the cation is within a defined cutoff.
        3. The distance between the cation and the aromatic ring system does not
        exceed a predefined distance cutoff.

    Args:
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.
        angle_cutoff (float): The maximum angle to consider for a cation pi contact.
        cache (bool): Determines whether computed results should be stored for later use to improve performance.

    Returns:
        NeighborPairs: The computed cation-pi contacts.

    """
    func = create_contact_function("cation_pi", angle_cutoff, cache)
    return func(ns)


def carbon_pi(ns: NeighborPairs, angle_cutoff: float = 30.0, cache: bool = False) -> NeighborPairs:
    """Handle the computation of carbon pi contacts in a molecular system.

    Carbon pi contacts are interactions involving the π system of carbon atoms in a molecule.
    This class, a derivative of the `ContactAnalysis` base class, overrides the `compute`
    method to provide functionality specifically for carbon pi contact computation.

    Carbon-pi contacts represent interactions between a carbon atom and an aromatic ring system.
    These contacts arise due to the partial positive charge on the carbon atom interacting with the electron-rich
    π-system of the aromatic ring.

    !!! tip "Definition"
        1. The interaction is between an aromatic ring and a carbon atom which is a weak hydrogen bond donor.
        2. The angle between the aromatic ring plane and the vector connecting the center of the aromatic ring
        and the carbon atom is within a specified cutoff.
        3. The distance between the carbon atom and the aromatic ring system does not
        exceed a predefined distance cutoff.

    Args:
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.
        angle_cutoff (float): The maximum angle to consider for a carbon pi contact.
        cache (bool): Determines whether computed results should be stored for later use to improve performance.

    Returns:
        NeighborPairs: The computed carbon-pi contacts.

    """
    func = create_contact_function("carbon_pi", angle_cutoff, cache)
    return func(ns)


def donor_pi(ns: NeighborPairs, angle_cutoff: float = 30.0, cache: bool = False) -> NeighborPairs:
    """Handle the computation of donor pi contacts in a molecular system.

    Donor pi contacts are interactions between electron donors and the π system of atoms in a molecule.
    This class, a derivative of the `ContactAnalysis` base class, overrides the `compute`
    method to provide functionality specifically for donor pi contact computation.

    Donor-pi contacts represent interactions between a hydrogen bond donor and an aromatic ring system.
    The electron-rich  aromatic ring can interact with the partial positive charge of a hydrogen atom
    that is involved in a polar bond, creating a stabilizing interaction.

    !!! tip "Definition"
        1. The interaction is between an aromatic ring and a hydrogen bond donor.
        2. The angle between the aromatic ring plane and the vector connecting the center of the aromatic ring
        and the hydrogen bond donor is within a predefined cutoff.
        3. The distance between the hydrogen bond donor and the aromatic ring system
            does not exceed a specified distance cutoff.

    Args:
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.
        angle_cutoff (float): The maximum angle to consider for a donor pi contact.
        cache (bool): Determines whether computed results should be stored for later use to improve performance.

    Returns:
        NeighborPairs: The computed donor-pi contacts.

    """
    func = create_contact_function("donor_pi", angle_cutoff, cache)
    return func(ns)


def sulphur_pi(ns: NeighborPairs, cache: bool = False) -> NeighborPairs:
    """Handle the computation of sulphur pi contacts in a molecular system.

    Sulphur pi contacts are interactions involving the π system of sulphur atoms in a molecule.
    This class, a derivative of the `ContactAnalysis` base class, overrides the `compute`
    method to provide functionality specifically for sulphur pi contact computation.

    Sulphur-pi interactions involve the interaction between an aromatic ring and a sulphur atom.
    Sulphur, particularly from methionine residues (MET), can interact with the electron cloud of an aromatic ring,
    contributing to the stability and specificity of biomolecular structures.

    !!! tip "Definition"
        1. The interaction is between an aromatic ring and a sulphur atom specifically in a methionine residue.
        2. The distance between the sulphur atom and the aromatic ring system
        does not exceed a predefined distance cutoff.

    Args:
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.
        cache (bool): Determines whether computed results should be stored for later use to improve performance.

    Returns:
        NeighborPairs: The computed sulphur-pi contacts.

    """
    func = create_contact_function("sulphur_pi", None, cache)
    return func(ns)


class AtomPlaneContacts:
    """Calculate and handles special atomic contacts within a molecular system, including
    carbon-pi, cation-pi, donor-pi, and sulphur-pi interactions. Each interaction type is computed
    as a method of this class.

    Aromatic rings and their interactions play a pivotal role in this analysis.

    !!! tip "Definition"
        1. The interaction is between an aromatic ring and a specific type of atom or group.
        2. The angle between the aromatic ring plane and the vector connecting the center of the
           aromatic ring and the specific atom or group is within a predefined cutoff (where applicable).
        3. The distance between the atom or group and the aromatic ring system does not exceed a
           predefined distance cutoff.

    The computation is based on the neighbors and angles calculated within the molecular system.

    Attributes:
        max_cutoff (float): The predefined maximum distance cutoff for atom-aromatic ring interactions.
        angles (Optional[NDArray[np.float32]]): Calculated angles between ring plane and atom vector.
        rings (Rings): Enumeration of rings in the molecular system.
        ap_contacts (_AtomPlaneContacts): Instance of `_AtomPlaneContacts` class.
        neighbors (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.

    ??? example "Example"
        ``` py
        luni = Luni(...)
        ns = luni.compute_neighbors()

        apc = AtomPlaneContacts(ns)
        dop = apc.donor_pi()    # for donor-pi contacts
        sup = apc.sulphur_pi()  # for sulphur-pi contacts
        cbp = apc.carbon_pi()   # for carbon-pi contacts
        cap = apc.cation_pi()   # for cation-pi contacts
        ```
    """

    max_cutoff = CONTACTS["aromatic"]["met_sulphur_aromatic_distance"]

    def __init__(self, ns: NeighborPairs):
        self.angles: Optional[NDArray[np.float32]] = None
        self.rings = enumerate_rings(ns.mol)
        self.ap_contacts = _AtomPlaneContacts()

        self._compute(ns, ns.luni.to("mda"))

    def _compute(self, ns: NeighborPairs, mda: AtomGroupType) -> None:
        pairs, distances = compute_neighbors(mda.atoms.positions, self.rings.centers)

        neighbors = ns.clone(pairs, distances)

        self.neighbors = neighbors - neighbors.type_filter("aromatic", partner=2)

        result = calc_ringnormal_pos_angle(self.neighbors, mda.universe.atoms, self.rings.centers, self.rings.normals)
        self.angles = result

    def donor_pi(self) -> NeighborPairs:
        """Handle the computation of donor pi contacts in a molecular system.

        Donor pi contacts are interactions between electron donors and the π system of atoms in a molecule.
        This class, a derivative of the `ContactAnalysis` base class, overrides the `compute`
        method to provide functionality specifically for donor pi contact computation.

        Donor-pi contacts represent interactions between a hydrogen bond donor and an aromatic ring system.
        The electron-rich  aromatic ring can interact with the partial positive charge of a hydrogen atom
        that is involved in a polar bond, creating a stabilizing interaction.

        !!! tip "Definition"
            1. The interaction is between an aromatic ring and a hydrogen bond donor.
            2. The angle between the aromatic ring plane and the vector connecting the center of the aromatic ring
            and the hydrogen bond donor is within a predefined cutoff.
            3. The distance between the hydrogen bond donor and the aromatic ring system
                does not exceed a specified distance cutoff.

        Returns:
            NeighborPairs: The computed donor-pi contacts.

        """
        assert self.angles is not None
        return self.ap_contacts.donor_pi(self.neighbors, self.angles)

    def sulphur_pi(self) -> NeighborPairs:
        """Handle the computation of sulphur pi contacts in a molecular system.

        Sulphur pi contacts are interactions involving the π system of sulphur atoms in a molecule.
        This class, a derivative of the `ContactAnalysis` base class, overrides the `compute`
        method to provide functionality specifically for sulphur pi contact computation.

        Sulphur-pi interactions involve the interaction between an aromatic ring and a sulphur atom.
        Sulphur, particularly from methionine residues (MET), can interact with the electron cloud of an aromatic ring,
        contributing to the stability and specificity of biomolecular structures.

        !!! tip "Definition"
            1. The interaction is between an aromatic ring and a sulphur atom specifically in a methionine residue.
            2. The distance between the sulphur atom and the aromatic ring system
            does not exceed a predefined distance cutoff.

        Returns:
            NeighborPairs: The computed sulphur-pi contacts.

        """
        assert self.angles is not None
        return self.ap_contacts.sulphur_pi(self.neighbors, self.angles)

    def carbon_pi(self) -> NeighborPairs:
        """Handle the computation of carbon pi contacts in a molecular system.

        Carbon pi contacts are interactions involving the π system of carbon atoms in a molecule.
        This class, a derivative of the `ContactAnalysis` base class, overrides the `compute`
        method to provide functionality specifically for carbon pi contact computation.

        Carbon-pi contacts represent interactions between a carbon atom and an aromatic ring system.
        These contacts arise due to the partial positive charge on the carbon atom interacting with the electron-rich
        π-system of the aromatic ring.

        !!! tip "Definition"
            1. The interaction is between an aromatic ring and a carbon atom which is a weak hydrogen bond donor.
            2. The angle between the aromatic ring plane and the vector connecting the center of the aromatic ring
            and the carbon atom is within a specified cutoff.
            3. The distance between the carbon atom and the aromatic ring system does not
            exceed a predefined distance cutoff.

        Returns:
            NeighborPairs: The computed carbon-pi contacts.

        """
        assert self.angles is not None
        return self.ap_contacts.carbon_pi(self.neighbors, self.angles)

    def cation_pi(self) -> NeighborPairs:
        """Handle the computation of cation pi contacts in a molecular system.

        Cation pi contacts are interactions between cations and the π system of atoms in a molecule.
        This class, a derivative of the `ContactAnalysis` base class, overrides the `compute`
        method to provide functionality specifically for cation pi contact computation.

        Cation-pi contacts involve interactions between a cation (a positively charged ion)
        and an electron-rich aromatic ring system. The aromatic system is capable of stabilizing the cation
        through its delocalized pi electrons.

        !!! tip "Definition"
            1. The interaction is between an aromatic ring and a positively ionizable atom.
            2. The angle between the aromatic ring plane and the vector connecting the center of the aromatic ring
            and the cation is within a defined cutoff.
            3. The distance between the cation and the aromatic ring system does not
            exceed a predefined distance cutoff.

        Returns:
            NeighborPairs: The computed cation-pi contacts.

        """
        assert self.angles is not None
        return self.ap_contacts.cation_pi(self.neighbors, self.angles)
