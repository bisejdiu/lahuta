from typing import Any, Callable, Dict, List, Tuple

import numpy as np
from joblib import Memory  # type: ignore
from MDAnalysis.lib import distances as mda_distances  # type: ignore
from numpy.typing import NDArray

from lahuta.config.defaults import CONTACTS
from lahuta.contacts.plane_plane import perceive_rings, vector_angle
from lahuta.core.neighbors import NeighborPairs
from lahuta.types.mda_commands import CappedDistance, DistanceType
from lahuta.types.mdanalysis import AtomGroupType

memory = Memory("cachedir", verbose=0)


DEFAULT_CONTACT_DISTS = {
    "cation_pi": CONTACTS["aromatic"]["atom_aromatic_distance"],
    "donor_pi": CONTACTS["aromatic"]["atom_aromatic_distance"],
    "sulphur_pi": CONTACTS["aromatic"]["met_sulphur_aromatic_distance"],
    "carbon_pi": CONTACTS["aromatic"]["atom_aromatic_distance"],
}


class _AtomPlaneContacts:
    @staticmethod
    def _donor_pi(
        ns: NeighborPairs, angles: NDArray[np.float_], angle_cutoff: float = 30.0
    ) -> NeighborPairs:
        """Compute the contacts between aromatic rings and the donor pi system."""
        distance = DEFAULT_CONTACT_DISTS["donor_pi"]
        return (
            ns.numeric_filter(angles, angle_cutoff)
            .distance_filter(distance)
            .type_filter("hbond_donor", partner=2)
        )

    @staticmethod
    def _sulphur_pi(ns: NeighborPairs, *_) -> NeighborPairs:
        """Compute the contacts between aromatic rings and the sulphur pi system."""
        distance = DEFAULT_CONTACT_DISTS["sulphur_pi"]
        indices = ns.partner2.select_atoms("resname MET and element S").indices  # type: ignore
        return ns.index_filter(indices, partner=2).distance_filter(distance)

    @staticmethod
    def _carbon_pi(
        ns: NeighborPairs, angles: NDArray[np.float_], angle_cutoff: float = 30.0
    ) -> NeighborPairs:
        """Compute the contacts between aromatic rings and the carbon pi system."""
        distance = DEFAULT_CONTACT_DISTS["carbon_pi"]
        return (
            ns.numeric_filter(angles, angle_cutoff)
            .index_filter(ns.partner2.select_atoms("element C").indices, partner=2)  # type: ignore
            .distance_filter(distance)
            .type_filter("weak_hbond_donor", partner=2)
        )

    @staticmethod
    def _cation_pi(
        ns: NeighborPairs, angles: NDArray[np.float_], angle_cutoff: float = 30.0
    ) -> NeighborPairs:
        """Compute the contacts between aromatic rings and the cation pi system."""
        distance = DEFAULT_CONTACT_DISTS["cation_pi"]
        return (
            ns.numeric_filter(angles, angle_cutoff)
            .distance_filter(distance)
            .type_filter("pos_ionisable", partner=2)
        )


@memory.cache  # type: ignore
def compute_neighbors(
    positions: NDArray[np.float_], rings: List[Dict[str, Any]]
) -> Tuple[NDArray[np.int_], NDArray[np.float_]]:
    max_cutoff = CONTACTS["aromatic"]["met_sulphur_aromatic_distance"]
    reference: NDArray[np.float_] = np.array([ring["center"] for ring in rings])

    wrapper: DistanceType = CappedDistance(mda_distances)
    pairs, distances = wrapper.capped_distance(
        reference, positions, max_cutoff, return_distances=True
    )
    return pairs, distances


@memory.cache  # type: ignore
def compute_angles(
    ns: NeighborPairs, uv_atoms: AtomGroupType, rings: List[Dict[str, Any]]
):
    ring_centers = np.array([ring["center"] for ring in rings])
    ring_normals = np.array([ring["normal"] for ring in rings])

    ring_centers = ring_centers[ns.pairs[:, 0]]
    ring_normals = ring_normals[ns.pairs[:, 0]]

    angles = vector_angle(
        ring_normals,
        ring_centers - uv_atoms[ns.pairs[:, 1]].positions,
    )
    return angles


def subtract_aromatic_neighbors(
    ns: NeighborPairs, pairs: NDArray[np.int_], distances: NDArray[np.float_]
):
    cloned_neighbors = ns.clone(pairs, distances)
    neighbors = cloned_neighbors - cloned_neighbors.type_filter("aromatic", partner=2)
    return neighbors


def compute_contacts(
    contact_fn: Callable[..., Any], angle_cutoff: float, use_cache: bool
):
    def wrapped(ns: NeighborPairs):
        mol = ns.luni.to("mol")
        mda: AtomGroupType = ns.luni.to("mda")
        rings = perceive_rings(mol)

        neighbors_fn = compute_neighbors if not use_cache else compute_neighbors.call
        result = neighbors_fn(mda.atoms.positions, rings)
        pairs, distances = result[0] if use_cache else result

        neighbors = subtract_aromatic_neighbors(ns, pairs, distances)

        angles_fn = compute_angles if not use_cache else compute_angles.call
        result = angles_fn(neighbors, mda.universe.atoms, rings)
        angles = result[0] if use_cache else result

        return contact_fn(neighbors, angles, angle_cutoff)

    return wrapped


def create_contact_function(contact_type, angle_cutoff, use_cache=True):
    func_name = f"_{contact_type}"
    contact_fn = getattr(_AtomPlaneContacts, func_name)
    return compute_contacts(contact_fn, angle_cutoff, use_cache)


# user-facing functions
def cation_pi(n, angle_cutoff: float = 30.0, cache=True):
    func = create_contact_function("cation_pi", angle_cutoff, cache)
    return func(n)


def carbon_pi(n, angle_cutoff: float = 30.0, cache=True):
    func = create_contact_function("carbon_pi", angle_cutoff, cache)
    return func(n)


def donor_pi(n, angle_cutoff: float = 30.0, cache=True):
    func = create_contact_function("donor_pi", angle_cutoff, cache)
    return func(n)


def sulphur_pi(n, cache=True):
    func = create_contact_function("sulphur_pi", None, cache)
    return func(n)


class AtomPlaneContacts:
    max_cutoff = CONTACTS["aromatic"]["met_sulphur_aromatic_distance"]

    def __init__(self, ns):
        self.angles = None
        self.rings = perceive_rings(ns.luni.to("mol"))
        self._ap_contacts = _AtomPlaneContacts()

        self._compute(ns, ns.luni.to("mda"))

    def _compute(self, ns, mda):
        result = compute_neighbors.call(mda.atoms.positions, self.rings)
        pairs, distances = result[0]

        neighbors = ns.clone(pairs, distances)

        self.neighbors = neighbors - neighbors.type_filter("aromatic", partner=2)

        result = compute_angles.call(self.neighbors, mda.universe.atoms, self.rings)
        self.angles = result[0]

    def donor_pi(self):
        return self._ap_contacts._donor_pi(self.neighbors, self.angles)

    def sulphur_pi(self):
        return self._ap_contacts._sulphur_pi(self.neighbors, self.angles)

    def carbon_pi(self):
        return self._ap_contacts._carbon_pi(self.neighbors, self.angles)

    def cation_pi(self):
        return self._ap_contacts._cation_pi(self.neighbors, self.angles)
