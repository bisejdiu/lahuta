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

Warning:
    While functionals maintain consistency in the contact computation API, 
    the class-based approach is more efficient as it avoids redundant computations.
    If speed is a priority, the class-based approach is recommended.

"""

import numpy as np
from MDAnalysis.lib import distances as mda_distances
from numpy.typing import NDArray

from lahuta.config.defaults import CONTACTS
from lahuta.core.neighbors import NeighborPairs
from lahuta.utils.array_utils import non_matching_indices, sorting_indices
from lahuta.utils.math import calc_vec_line_angles
from lahuta.utils.ob import enumerate_rings

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

def compute_neighbors(
    positions: NDArray[np.float32], reference: NDArray[np.float32]
) -> tuple[NDArray[np.int32], NDArray[np.float32]]:
    """Compute the neighbors between the reference and positions."""
    max_cutoff = CONTACTS["aromatic"]["met_sulphur_aromatic_distance"]

    pairs, distances = mda_distances.capped_distance(reference, positions, max_cutoff, return_distances=True)
    return pairs, distances

def calc_ringnormal_pos_angle(
    positions: NDArray[np.float32], ring_centers: NDArray[np.float32], ring_normals: NDArray[np.float32]
) -> NDArray[np.float32]:
    """Calculate the angle between the ring normal and the vector connecting the ring center and the atom."""
    return calc_vec_line_angles(
        ring_normals,
        ring_centers - positions,
    )


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
        angles (Optional[NDArray[np.float32]]): Calculated angles between ring plane and atom vector.
        rings (Rings): Enumeration of rings in the molecular system.
        ns (NeighborPairs): The object encapsulating pairs of neighboring atoms in the system.
        
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

    def __init__(self, ns: NeighborPairs):
        mda = ns.luni.to("mda")

        self.rings = enumerate_rings(ns.luni.to("mol"))
        self.ring_first_idxs = mda.universe.atoms.indices[self.rings.first_atom_idx]

        pairs, distances = compute_neighbors(mda.atoms.positions, self.rings.centers)
        mapped_pairs = np.array([self.ring_first_idxs[pairs[:, 0]], pairs[:, 1]]).T

        nn = NeighborPairs(ns.luni)
        nn.set_neighbors(mapped_pairs, distances, sort=False)

        npairs = non_matching_indices(mapped_pairs, nn.type_filter("aromatic", partner=2).pairs)
        pairs = pairs[npairs]
        nn.set_neighbors(nn.pairs[npairs], nn.distances[npairs], sort=False)

        indices = sorting_indices(nn.pairs)
        nn.set_neighbors(nn.pairs[indices], nn.distances[indices], sort=False)
        pairs = pairs[indices]

        self.ns = nn

        self.angles = calc_ringnormal_pos_angle(
            mda.universe.atoms[pairs[:, 1]].positions, 
            self.rings.centers[pairs[:, 0]], 
            self.rings.normals[pairs[:, 0]]
        )


        # mapped_neighbors = ns.new(mapped_pairs, distances)
        # neighbors = ns.new(pairs, distances)

        # self.ns_ = neighbors - neighbors.type_filter("aromatic", partner=2)
        # self.ns = mapped_neighbors - mapped_neighbors.type_filter("aromatic", partner=2)

        # self.angles = calc_ringnormal_pos_angle2(
        #     mda.universe.atoms[self.ns_.pairs[:, 1]].positions, 
        #     self.rings.centers[self.ns_.pairs[:, 0]], 
        #     self.rings.normals[self.ns_.pairs[:, 0]]
        # )

    def donor_pi(self, angle_cutoff: float = 30.0) -> NeighborPairs:
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

        Returns:
            NeighborPairs: The computed donor-pi contacts.
        """
        distance = DEFAULT_CONTACT_DISTS["donor_pi"]
        return (
            self.ns.numeric_filter(self.angles, angle_cutoff)
            .distance_filter(distance)
            .type_filter("hbond_donor", partner=2)
        )
    
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
        distance = DEFAULT_CONTACT_DISTS["sulphur_pi"]
        indices = self.ns.partner2.select_atoms("resname MET and element S").indices
        return self.ns.index_filter(indices, partner=2).distance_filter(distance)
    
    def carbon_pi(self, angle_cutoff: float = 30.0) -> NeighborPairs:
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

        Returns:
            NeighborPairs: The computed carbon-pi contacts.

        """
        distance = DEFAULT_CONTACT_DISTS["carbon_pi"]
        return (
            self.ns.numeric_filter(self.angles, angle_cutoff)
            .index_filter(self.ns.partner2.select_atoms("element C").indices, partner=2)
            .distance_filter(distance)
            .type_filter("weak_hbond_donor", partner=2)
        )
    
    def cation_pi(self, angle_cutoff: float = 30.0) -> NeighborPairs:
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

        Returns:
            NeighborPairs: The computed cation-pi contacts.

        """
        distance = DEFAULT_CONTACT_DISTS["cation_pi"]
        return (
            self.ns.numeric_filter(self.angles, angle_cutoff)
            .distance_filter(distance)
            .type_filter("pos_ionisable", partner=2)
        )
    



# user-facing functions
def cation_pi(ns: NeighborPairs, angle_cutoff: float = 30.0) -> NeighborPairs:
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
    apc = AtomPlaneContacts(ns)
    return apc.cation_pi(angle_cutoff)

def carbon_pi(ns: NeighborPairs, angle_cutoff: float = 30.0) -> NeighborPairs:
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
    apc = AtomPlaneContacts(ns)
    return apc.carbon_pi(angle_cutoff)

def donor_pi(ns: NeighborPairs, angle_cutoff: float = 30.0) -> NeighborPairs:
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
    apc = AtomPlaneContacts(ns)
    return apc.donor_pi(angle_cutoff)

def sulphur_pi(ns: NeighborPairs) -> NeighborPairs:
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
    apc = AtomPlaneContacts(ns)
    return apc.sulphur_pi()
