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

    if reference.shape[0] == 0 or positions.shape[0] == 0:
        return np.empty((0, 2), dtype=np.int32), np.empty(0, dtype=np.float32)
    pairs, distances = mda_distances.capped_distance(reference, positions, max_cutoff, return_distances=True)
    print("PAIRS size: ", pairs.shape)
    return pairs, distances**2


def calc_ringnormal_pos_angle(
    positions: NDArray[np.float32], ring_centers: NDArray[np.float32], ring_normals: NDArray[np.float32]
) -> NDArray[np.float32]:
    """Calculate the angle between the ring normal and the vector connecting the ring center and the atom."""
    # for i in range(10):
    #     print("positions: ", positions[i], "ring_centers: ", ring_centers[i], "ring_normals: ", ring_normals[i])
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

    def __init__(self, ns: NeighborPairs, allow_metal_participation: bool = False) -> None:
        # allow_metal_participation: A boolean option to include or exclude metal atoms in pi-system interactions.
        mda = ns.luni.to("mda")

        self.rings = enumerate_rings(ns.luni.to("mol"))
        print("____" * 10)
        # for ix, ring in enumerate(self.rings.rings):
        #     print("ring: ", ring, self.rings._atoms)
        # print("ring: ", self.rings._atoms)
        _atoms = self.rings._atoms
        for x in _atoms:
            if 1197 in x:  # or 1517 in x:
                xx = " ".join([str(y) for y in x])
                print("x: ", xx)

        print("number of rings: ", self.rings.centers.shape[0])
        _rings = ns.luni._file_loader.luni.get_rings()
        _first_atom_idx = []
        for _ring in _rings.rings:
            atom_ids = _ring.atom_ids
            atom_ids.sort()
            _first_atom_idx.append(atom_ids[0])

        self.ring_first_idxs = mda.universe.atoms.indices[self.rings.first_atom_idx]
        _first_atom_idx = np.array(_first_atom_idx)
        # print("Diff: ", np.all(self.ring_first_idxs == self.rings.first_atom_idx))
        # print("first_atom_idx: ", sorted(self.ring_first_idxs), type(self.ring_first_idxs))
        # print("_first_atom_idx: ", sorted(_first_atom_idx), type(_first_atom_idx))

        # no rings means no contacts
        if self.rings.centers.shape[0] == 0:
            self.ns = NeighborPairs(ns.luni)
            self.angles = np.empty(0, dtype=np.float32)
            return

        # FIX: A big issue is that the registered distances are to ring centers, but stored as if
        # they were atom-atom distances. This is very confusing!
        pairs, distances = compute_neighbors(mda.positions, self.rings.centers)
        _pairs, _distances = compute_neighbors(mda.positions, _rings.centers)

        # print("shape: ", mda.universe.atoms.positions.shape, self.rings.centers.shape, self.rings.normals.shape)
        # angs = calc_ringnormal_pos_angle(mda.universe.atoms.positions[:10], self.rings.centers, self.rings.normals)
        # print("angles: ", angs)

        # key: kurrent index, value: use this index
        # ring_order_indices = {
        #     2: 0,
        #     6: 1,
        #     3: 2,
        #     7: 3,
        #     8: 4,
        #     9: 5,
        #     5: 6,
        #     0: 7,
        #     1: 8,
        #     4: 9,
        # }
        # reversed_ring_order_indices = {v: k for k, v in ring_order_indices.items()}

        # print("computing 1 angle for each ring: ")
        # for ix, ring in enumerate(self.rings.rings):
        #     # print("ring: ", ring)
        #     pos = mda.universe.atoms[ix].position.reshape(1, 3)
        #     # print("pos: ", pos)
        #     use_index = ring_order_indices[ix]
        #     use_ring_center = self.rings.centers[use_index].reshape(1, 3)
        #     use_ring_normal = self.rings.normals[use_index].reshape(1, 3)
        #     print("center: ", self.rings.centers[ix])
        #     # print("normal: ", self.rings.normals[ix])
        #     ang = calc_ringnormal_pos_angle(
        #         pos,
        #         use_ring_center,
        #         use_ring_normal,
        #         # pos, self.rings.centers[ix].reshape(1, 3), self.rings.normals[ix].reshape(1, 3)
        #     )
        #     # print("angle: ", ang)

        # for k, v in reversed_ring_order_indices.items():
        #     for ix, ring in enumerate(self.rings.rings):
        #         if ix == v:
        #             pos = mda.universe.atoms[k].position.reshape(1, 3)
        #             # use_index = ring_order_indices[ix]
        #             use_ring_center = self.rings.centers[ix].reshape(1, 3)
        #             use_ring_normal = self.rings.normals[ix].reshape(1, 3)
        #             print("center: ", use_ring_center)
        #             print("normal: ", use_ring_normal)
        #             print("pos: ", pos)
        #             ang = calc_ringnormal_pos_angle(
        #                 pos,
        #                 use_ring_center,
        #                 use_ring_normal,
        #                 # pos, self.rings.centers[ix].reshape(1, 3), self.rings.normals[ix].reshape(1, 3)
        #             )
        #             print("angle: ", ang)

        NN = NeighborPairs(ns.luni)
        NN.set_neighbors(_pairs, _distances, sort=False)

        mapped_pairs = np.array([self.ring_first_idxs[pairs[:, 0]], pairs[:, 1]]).T
        _mapped_pairs = np.array([_first_atom_idx[_pairs[:, 0]], _pairs[:, 1]]).T

        # print("mapped_pairs: ", mapped_pairs)
        # print("_mapped_pairs: ", _mapped_pairs)

        nn = NeighborPairs(ns.luni)
        nn.set_neighbors(mapped_pairs, distances, sort=False)

        _nn = NeighborPairs(ns.luni)
        _nn.set_neighbors(_mapped_pairs, _distances, sort=False)

        # NOTE: negative selection for non aromatic atoms in the second partner
        # that's because this class handles atom-plane contacts, and we already have planes
        # in the rings system.
        npairs = non_matching_indices(mapped_pairs, nn.type_filter("aromatic", partner=2).pairs)
        pairs = pairs[npairs]

        _npairs = non_matching_indices(_mapped_pairs, _nn.type_filter("aromatic", partner=2).pairs)
        _pairs = _pairs[_npairs]

        nn.set_neighbors(nn.pairs[npairs], nn.distances[npairs], sort=False)
        _nn.set_neighbors(_nn.pairs[_npairs], _nn.distances[_npairs], sort=False)
        # for x in _nn:
        #     print("->", x.labels, np.sqrt(x.distances))

        from lahuta.config.atoms import METALS

        metal_indices = _nn.atoms.select_atoms("element " + " ".join(METALS)).indices
        self.metal_indices = metal_indices

        # neg_nn = _nn.index_filter(metal_indices, partner=1)
        # _nn = _nn - neg_nn

        # _nn.set_neighbors(tmp.pairs, tmp.distances, sort=False)
        # print("p shape: ", _nn.pairs.shape)
        # neg_nn = _nn.index_filter(metal_indices, partner=2)
        # print("neg_nn: ", neg_nn)
        # _nn = _nn - neg_nn

        indices = sorting_indices(nn.pairs)
        nn.set_neighbors(nn.pairs[indices], nn.distances[indices], sort=False)
        pairs = pairs[indices]

        _indices = sorting_indices(_nn.pairs)
        _nn.set_neighbors(_nn.pairs[_indices], _nn.distances[_indices], sort=False)
        # for x in _nn:
        #     print("->", x.pairs, x.distances)
        # for ix, x in enumerate(_nn.labels):
        #     print("->", x, np.sqrt(_nn.distances[ix]))
        _pairs = _pairs[_indices]

        # print("final nn: ", nn)

        self.ns = nn
        self._ns = _nn

        self.angles = calc_ringnormal_pos_angle(
            mda.universe.atoms[pairs[:, 1]].positions, self.rings.centers[pairs[:, 0]], self.rings.normals[pairs[:, 0]]
        )
        # print("angles: ", self.angles)
        # print("angles: ", self.angles.shape, np.unique(self.angles).shape)
        _angles = calc_ringnormal_pos_angle(
            mda.universe.atoms[_pairs[:, 1]].positions, _rings.centers[_pairs[:, 0]], _rings.norm1[_pairs[:, 0]]
        )
        self._angles = _angles
        # print("angles: ", self.angles.shape)
        # print("_angles: ", _angles.shape)

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
        # print("distance", distance)
        distance *= distance
        # print("distance", distance)
        # print("x", self.ns.distances)
        return (
            self._ns.numeric_filter(self._angles, angle_cutoff)
            .distance_filter(distance)
            .type_filter("hbond_donor", partner=2)
        )
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
        distance *= distance
        indices = self._ns.partner2.select_atoms("resname MET and element S").indices
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
        print("distance: ", distance)
        distance *= distance
        print("distance: ", distance)
        return (
            self._ns.numeric_filter(self._angles, angle_cutoff)
            .index_filter(self.ns.partner2.select_atoms("element C").indices, partner=2)
            .distance_filter(distance)
            .type_filter("weak_hbond_donor", partner=2)
            # .index_filter(self.metal_indices, partner=1, reverse=True)
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
        distance *= distance
        return (
            self._ns.numeric_filter(self._angles, angle_cutoff)
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
    # print("apc: ", apc.ns.distances)
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
