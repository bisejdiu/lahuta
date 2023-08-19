"""Provides both class and function-level APIs to compute plane-plane contacts
between two ring-containing residues or ligands. Each interaction is labeled according 
to the relative orientations of the planes.

Types of Plane-Plane Contact Orientations:
    FF (Face-Face): Represents two planes interacting in a face-to-face orientation, where 
        their flat surfaces are parallel and directly aligned. This typically leads to 
        maximum pi-pi stacking.
    OF (Oblique Face): Denotes a skewed face-to-face orientation. One plane is tilted 
        relative to the other, resulting in partial overlap of their surfaces.
    EE (Edge-Edge): Describes an interaction where the edges of two planes align parallel, 
        but their surfaces do not overlap, limiting pi-pi interactions.
    FT (Face-T-stack): Indicates an interaction where the face of one plane (the base of the 'T') 
        aligns with the edge of another plane (the stem of the 'T'). The second plane is 
        oriented in a T-stacked manner.
    OT (Oblique-T-stack): Represents an interaction where an obliquely oriented face of one 
        plane (the skewed base of the 'T') aligns with the edge of another plane (the stem of 
        the 'T'). The second plane is in a T-stacked orientation.
    ET (Edge-T-stack): Describes an interaction where the edge of one plane aligns with the 
        edge of another plane. Both are part of the stem in a T-stacked orientation.
    FE (Face-Edge): Depicts an interaction where the flat surface (face) of one plane aligns 
        with the edge of another plane, leading to partial pi-pi interactions.
    OE (Oblique-Edge): Represents an interaction where an obliquely oriented face of one plane 
        aligns with the edge of another plane, leading to partial and skewed pi-pi interactions.
    EF (Edge-Face): Indicates an interaction where the edge of one plane aligns with the face of 
        another plane, creating limited pi-pi interactions.

Classes:
    PlanePlaneContacts: A class to compute plane-plane contacts.

Functions:
    plane_plane_neighbors: Function to compute plane-plane contacts.

Usage:
    luni = Luni(...)
    ns = luni.compute_neighbors()
    plane_plane = PlanePlaneContacts(ns)
    result = plane_plane.results
"""

from typing import Any, Dict, List, Tuple

import numpy as np
from numpy.typing import NDArray

from lahuta.config.defaults import CONTACTS
from lahuta.core.neighbors import NeighborPairs
from lahuta.utils.array_utils import sorting_indices
from lahuta.utils.math import calc_vec_angle
from lahuta.utils.ob import enumerate_rings


class _PlanePlaneContacts:
    def __init__(self, ns: NeighborPairs) -> None:
        self.ns = ns
        self.rings = enumerate_rings(self.ns.mol)
        self.centroid_distance = CONTACTS["aromatic"]["centroid_distance"]
        self._annotations: Dict[str, NDArray[Any]] = {}

        self._pair_ids: NDArray[np.int32] = np.array([])
        self.distances: NDArray[np.float32] = np.array([])

    def compute(self) -> None:
        """Compute plane-plane contacts based on the neighbor pairs and the set centroid distance."""
        centers = self.rings.centers
        normals = self.rings.normals

        pair_ids = self._gen_combinations()

        # compute the pairwise distances
        pair_dists = np.linalg.norm(centers[pair_ids[:, 0]] - centers[pair_ids[:, 1]], axis=1)

        # apply the distance mask
        valid_distance_mask = pair_dists <= self.centroid_distance
        pair_ids = pair_ids[valid_distance_mask]
        pair_dists = pair_dists[valid_distance_mask]

        first_pair_ids, second_pair_ids = pair_ids[:, 0], pair_ids[:, 1]
        pair_diffs = centers[first_pair_ids] - centers[second_pair_ids]
        normal_angles = calc_vec_angle(normals[first_pair_ids], normals[second_pair_ids])
        theta_angles = calc_vec_angle(normals[first_pair_ids], pair_diffs)

        int_types = assign_pp_contact_type(normal_angles, theta_angles)

        # check for int_type is not None
        # TODO @bisejdiu: Shouldn't this be None, without quotes?
        mask = int_types != "None"

        self._pair_ids = pair_ids[mask]
        self.distances = pair_dists[mask]

        self._make_annotations(theta_angles[mask], normal_angles[mask], int_types[mask])

    def _make_annotations(
        self,
        theta_angles: NDArray[np.float32],
        normal_angles: NDArray[np.float32],
        int_types: NDArray[np.str_],
    ) -> None:
        ring_atoms = self.rings.atoms[self._pair_ids]
        self._annotations["theta_angles"] = theta_angles
        self._annotations["normal_angles"] = normal_angles
        self._annotations["ring1_atoms"] = ring_atoms[:, 0]
        self._annotations["ring2_atoms"] = ring_atoms[:, 1]
        self._annotations["contact_labels"] = int_types

    def _get_pairs_distances(self) -> Tuple[NDArray[np.int32], NDArray[np.float32]]:
        ring_atom_indices = self.rings.first_atom_idx[self._pair_ids]
        first_ring_indices, second_ring_indices = (
            ring_atom_indices[:, 0],
            ring_atom_indices[:, 1],
        )

        pairs = np.array(list(zip(first_ring_indices, second_ring_indices)))

        return pairs, self.distances

    def _sort_inputs(self) -> Tuple[NDArray[np.int32], NDArray[np.float32]]:
        pairs, distances = self._get_pairs_distances()
        indices_arr = sorting_indices(pairs)
        indices: List[int] = indices_arr.tolist()
        pairs, distances = NeighborPairs.sort_inputs(pairs, self.distances)

        self._pair_ids = self._pair_ids[indices]
        sorted_annotations = {k: v[indices] for k, v in self._annotations.items()}
        self._annotations = sorted_annotations

        return pairs, distances

    def get_neighbors(self) -> NeighborPairs:
        """Return the plane-plane contacts as a NeighborPairs object."""
        pairs, distances = self._sort_inputs()
        ns = self.ns.clone(pairs, distances)
        ns.annotations = self._annotations
        return ns

    def _gen_combinations(self, use_itertools: bool = False) -> NDArray[np.int32]:
        """Generate all combinations of pairs of indices in the form (i, j) where i < j."""
        if use_itertools:
            from itertools import combinations

            return np.array(list(combinations(range(len(self.rings)), 2)))

        # scales better with the number of rings
        return np.column_stack(np.triu_indices(len(self.rings), k=1))


def plane_plane_neighbors(ns: NeighborPairs) -> NeighborPairs:
    """Handle the computation of plane-plane contacts in a molecular system.

    This class provides both class and function-level APIs to compute plane-plane contacts
    between two ring-containing residues or ligands. Each interaction is labeled according
    to the relative orientations of the planes.

    !!! tip "Definitions"
        **FF** _(Face-Face)_: Represents two planes interacting in a face-to-face orientation, where
            their flat surfaces are parallel and directly aligned. This typically leads to
            maximum pi-pi stacking.

        **OF** _(Oblique Face)_: Denotes a skewed face-to-face orientation. One plane is tilted
            relative to the other, resulting in partial overlap of their surfaces.

        **EE** _(Edge-Edge)_: Describes an interaction where the edges of two planes align parallel,
            but their surfaces do not overlap, limiting pi-pi interactions.

        **FT** _(Face-T-stack)_: Indicates an interaction where the face of one plane (the base of the 'T')
            aligns with the edge of another plane (the stem of the 'T'). The second plane is
            oriented in a T-stacked manner.

        **OT** _(Oblique-T-stack)_: Represents an interaction where an obliquely oriented face of one
            plane (the skewed base of the 'T') aligns with the edge of another plane (the stem of
            the 'T'). The second plane is in a T-stacked orientation.

        **ET** _(Edge-T-stack)_: Describes an interaction where the edge of one plane aligns with the
            edge of another plane. Both are part of the stem in a T-stacked orientation.

        **FE** _(Face-Edge)_: Depicts an interaction where the flat surface (face) of one plane aligns
            with the edge of another plane, leading to partial pi-pi interactions.

        **OE** _(Oblique-Edge)_: Represents an interaction where an obliquely oriented face of one plane
            aligns with the edge of another plane, leading to partial and skewed pi-pi interactions.

        **EF** _(Edge-Face)_: Indicates an interaction where the edge of one plane aligns with the face of
            another plane, creating limited pi-pi interactions.

    Args:
        ns (NeighborPairs): A NeighborPairs object containing the atom neighbor relationships in the system.

    Returns:
        NeighborPairs: A NeighborPairs object containing only plane-plane contacts.

    """
    plane_plane = _PlanePlaneContacts(ns)
    plane_plane.compute()
    return plane_plane.get_neighbors()


class PlanePlaneContacts:
    """Handle the computation of plane-plane contacts in a molecular system.

    This class provides both class and function-level APIs to compute plane-plane contacts
    between two ring-containing residues or ligands. Each interaction is labeled according
    to the relative orientations of the planes.

    !!! tip "Definitions"
        **FF** _(Face-Face)_: Represents two planes interacting in a face-to-face orientation, where
            their flat surfaces are parallel and directly aligned. This typically leads to
            maximum pi-pi stacking.

        **OF** _(Oblique Face)_: Denotes a skewed face-to-face orientation. One plane is tilted
            relative to the other, resulting in partial overlap of their surfaces.

        **EE** _(Edge-Edge)_: Describes an interaction where the edges of two planes align parallel,
            but their surfaces do not overlap, limiting pi-pi interactions.

        **FT** _(Face-T-stack)_: Indicates an interaction where the face of one plane (the base of the 'T')
            aligns with the edge of another plane (the stem of the 'T'). The second plane is
            oriented in a T-stacked manner.

        **OT** _(Oblique-T-stack)_: Represents an interaction where an obliquely oriented face of one
            plane (the skewed base of the 'T') aligns with the edge of another plane (the stem of
            the 'T'). The second plane is in a T-stacked orientation.

        **ET** _(Edge-T-stack)_: Describes an interaction where the edge of one plane aligns with the
            edge of another plane. Both are part of the stem in a T-stacked orientation.

        **FE** _(Face-Edge)_: Depicts an interaction where the flat surface (face) of one plane aligns
            with the edge of another plane, leading to partial pi-pi interactions.

        **OE** _(Oblique-Edge)_: Represents an interaction where an obliquely oriented face of one plane
            aligns with the edge of another plane, leading to partial and skewed pi-pi interactions.

        **EF** _(Edge-Face)_: Indicates an interaction where the edge of one plane aligns with the face of
            another plane, creating limited pi-pi interactions.

    Args:
        ns (NeighborPairs): A NeighborPairs object containing the atom neighbor relationships in the system.

    Attributes:
        centroid_distance (float): The maximum distance to consider for the centroid of a plane-plane contact.
            See `lahuta.config.defaults.CONTACTS` for default values.
        results (NeighborPairs): NeighborPairs containing the pairs of atoms found to be forming plane-plane contacts.

    ??? example "Example"
        ``` py
        luni = Luni(...)
        ns = luni.compute_neighbors()

        plane_plane = PlanePlaneContacts(ns)
        result = plane_plane.results
        ```
    """

    centroid_distance = CONTACTS["aromatic"]["centroid_distance"]

    def __init__(self, ns: NeighborPairs) -> None:
        self.ns = ns
        self.results = self.compute()

    def compute(self) -> NeighborPairs:
        """Compute plane-plane contacts based on the neighbor pairs and the set centroid distance.

        This method initializes a `_PlanePlaneContacts` object with the neighbor pairs object (`ns`),
        sets its `centroid_distance` to match the class's `centroid_distance`, calls its `compute`
        method, and finally calls its `get_neighbors` method to obtain the plane-plane contacts.

        Returns
        -------
        NeighborPairs
            The object containing the pairs of atoms that are considered as plane-plane contacts
            based on the set centroid_distance criteria.
        """
        pp = _PlanePlaneContacts(self.ns)
        pp.centroid_distance = self.centroid_distance
        pp.compute()

        return pp.get_neighbors()


def assign_pp_contact_type(normal_angle: NDArray[np.float32], theta: NDArray[np.float32]) -> NDArray[np.str_]:
    """Assign a contact type based on the normal angle and theta angle.

    This function assigns a contact type based on the normal angle and theta angle between two planes.
    The normal angle is the angle between the normals of the two planes, while the theta angle is the
    angle between the normal of the first plane and the vector connecting the centroids of the two planes.

    Args:
        normal_angle (NDArray[np.float32]): The angle between the normals of the two planes.
        theta (NDArray[np.float32]): The angle between the normal of the first plane and the vector
            connecting the centroids of the two planes.

    Returns:
        NDArray[np.str_]: An array of strings representing the contact type for each pair of planes.

    """
    types = np.array(["FF", "OF", "EE", "FT", "OT", "ET", "FE", "OE", "EF"], dtype=np.str_)
    ranges = np.array(
        [
            (0, 30, 0, 30),
            (0, 30, 30, 60),
            (0, 30, 60, 90),
            (30, 60, 0, 30),
            (30, 60, 30, 60),
            (30, 60, 60, 90),
            (60, 90, 0, 30),
            (60, 90, 30, 60),
            (60, 90, 60, 90),
        ],
        dtype=np.int32,
    )

    conditions = np.logical_and(
        np.logical_and(ranges[:, 0] <= normal_angle[:, None], normal_angle[:, None] <= ranges[:, 1]),
        np.logical_and(ranges[:, 2] <= theta[:, None], theta[:, None] <= ranges[:, 3]),
    )

    indices: NDArray[np.int32] = np.argmax(conditions, axis=1)
    return types[indices]
