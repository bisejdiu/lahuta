"""
Module: plane_plane.py

This module provides both class and function-level APIs to compute plane-plane contacts 
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
    universe = Universe(...)
    ns = universe.compute_neighbors()
    plane_plane = PlanePlaneContacts(ns)
    result = plane_plane.compute().results
"""

from typing import Any, Dict, List, Tuple

import numpy as np
from numpy.typing import NDArray
from openbabel import openbabel as ob

from lahuta.config.defaults import CONTACTS
from lahuta.core.neighbors import NeighborPairs
from lahuta.lahuta_types.openbabel import MolType
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
        mask = int_types != "None"  # FIXME: Shouldn't this be None, without quotes?

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
        indices_arr = NeighborPairs.get_sorting_index(pairs)
        indices: List[int] = indices_arr.tolist()
        pairs, distances = NeighborPairs.sort_inputs(pairs, self.distances)

        self._pair_ids = self._pair_ids[indices]
        sorted_annotations = {k: v[indices] for k, v in self._annotations.items()}
        self._annotations = sorted_annotations

        return pairs, distances

    def get_neighbors(self) -> NeighborPairs:
        pairs, distances = self._sort_inputs()
        ns = self.ns.clone(pairs, distances)
        ns.annotations = self._annotations
        return ns

    def _gen_combinations(self, use_itertools: bool = False) -> NDArray[np.int32]:
        """Generate all combinations of pairs of indices in the form (i, j) where i < j"""
        if use_itertools:
            # pylint: disable=import-outside-toplevel
            from itertools import combinations  # type: ignore

            return np.array(list(combinations(range(len(self.rings)), 2)))
        else:
            # scales better with the number of rings
            return np.column_stack(np.triu_indices(len(self.rings), k=1))


def plane_plane_neighbors(ns: NeighborPairs) -> NeighborPairs:
    pp = _PlanePlaneContacts(ns)
    pp.compute()
    return pp.get_neighbors()


class PlanePlaneContacts:
    """Handles the computation of plane-plane contacts in a molecular system.

    This class provides both class and function-level APIs to compute plane-plane contacts
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

    Attributes
    ----------
    ns : NeighborPairs
        The object encapsulating pairs of neighboring atoms in the system.
        This object is usually generated by calling the `compute_neighbors` method
        on the Lahuta Universe object.
    centroid_distance : float
        The maximum distance to consider for the centroid of a plane-plane contact.
        This value is retrieved from the 'centroid_distance' entry of the 'aromatic'
        section in the global CONTACTS dictionary.
    results : NeighborPairs
        The object containing the pairs of atoms that are considered as plane-plane contacts
        based on the set centroid_distance criteria. It is computed when the class is instantiated.

    Methods
    -------
    compute() -> NeighborPairs:
        Computes and returns the plane-plane contacts, utilizing the `_PlanePlaneContacts` class and
        setting its `centroid_distance` to match the class's `centroid_distance`.
    """

    centroid_distance = CONTACTS["aromatic"]["centroid_distance"]

    def __init__(self, ns: NeighborPairs) -> None:
        self.ns = ns
        self.results = self.compute()

    def compute(self) -> NeighborPairs:
        """Computes plane-plane contacts based on the neighbor pairs and the set centroid distance.

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


# def vector_angle_1d(v1: NDArray[np.float32], v2: NDArray[np.float32]) -> NDArray[np.float32]:
#     v1_u = v1 / np.linalg.norm(v1, axis=-1, keepdims=True)
#     v2_u = v2 / np.linalg.norm(v2, axis=-1, keepdims=True)
#     dot = np.einsum("ij,ij->i", v1_u, v2_u)  # type: ignore
#     angle: NDArray[np.float32] = np.arccos(np.clip(dot, -1.0, 1.0))
#     angle = np.sign(dot) * angle
#     angle_deg: NDArray[np.float32] = np.degrees(angle)
#     angle_deg[angle_deg < 0] = 180 + angle_deg[angle_deg < 0]
#     return angle_deg


def assign_pp_contact_type(normal_angle: NDArray[np.float32], theta: NDArray[np.float32]) -> NDArray[np.str_]:
    """Assigns a contact type based on the normal angle and theta angle."""
    types = np.array(["FF", "OF", "EE", "FT", "OT", "ET", "FE", "OE", "EF"])
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
        ]
    )

    conditions = np.logical_and(
        np.logical_and(ranges[:, 0] <= normal_angle[:, None], normal_angle[:, None] <= ranges[:, 1]),
        np.logical_and(ranges[:, 2] <= theta[:, None], theta[:, None] <= ranges[:, 3]),
    )

    indices = np.argmax(conditions, axis=1)
    return types[indices]  # type: ignore


# def vector_angle(v1: NDArray[np.float_], v2: NDArray[np.float_]) -> NDArray[np.float_]:
#     """Calculate the angle between two vectors.

#     Parameters
#     ----------
#     v1 : np.ndarray
#         First vector. It must be already normalized.
#     v2 : np.ndarray
#         Second vector.

#     Returns
#     -------
#     float
#         Angle between vectors in degrees.
#     """

#     dot = np.einsum("ij,ij->i", v1, v2 / np.linalg.norm(v2, axis=1)[..., np.newaxis])  # type: ignore
#     angle = np.sign(dot) * np.arccos(dot)

#     angle = np.where(angle < 0, angle + np.pi, angle)

#     return np.degrees(angle)  # type: ignore


# FIXME: atom_plane needs to use the Rings class, rather than `perceive_rings`
def perceive_rings(mol: MolType) -> List[Dict[str, Any]]:
    aromatics: List[Dict[str, Any]] = []

    for _, ob_ring in enumerate(mol.GetSSSR()):
        if not ob_ring.IsAromatic():
            continue

        center = ob.vector3()
        normal = ob.vector3()

        ob_ring.findCenterAndNormal(center, normal, ob.vector3())

        center = np.array([center.GetX(), center.GetY(), center.GetZ()])  # type: ignore
        normal = np.array([normal.GetX(), normal.GetY(), normal.GetZ()])  # type: ignore
        atoms = sorted([atom for atom in ob_ring._path])  # type: ignore
        aromatics.append(
            {
                "atoms": atoms,
                "center": center,
                "normal": normal,
            }
        )

    return aromatics


# def process_ring_information(indices, rings):
#     """Process ring information for a given set of indices."""
#     ring_atoms = {ix: ring["atoms"] for ix, ring in enumerate(rings) if ix in indices}
#     ring_atoms = [ring_atoms[ix] for ix in indices]
#     first_ring_atoms = [ring[0] for ring in ring_atoms]

#     return first_ring_atoms, ring_atoms


# class APDataFrameFactory:
#     """A class for storing and manipulating contact data."""

#     def __init__(
#         self,
#         apcontacts: Any,
#         rings: Any,
#         label: str = "ring",
#         df_format: Literal["print", "compact", "expanded"] = "print",
#     ):
#         """A class for storing and manipulating contact data."""

#         nag = apcontacts._atoms[apcontacts.pairs]

#         first_rings, ring_atoms = process_ring_information(nag[:, 0].indices, rings)

#         col1, col2 = nag.universe.atoms[first_rings], nag[:, 1]

#         self.methods = ["resids", "resnames", "names", "indices"]
#         self.format = df_format

#         data = {}
#         for method in self.methods:
#             data[f"residue1_{method}"] = getattr(col1, method)
#             data[f"residue2_{method}"] = getattr(col2, method)

#         data["ring_atoms"] = ring_atoms
#         data["distances"] = apcontacts.distances
#         data["labels"] = [label] * len(apcontacts.distances)

#         self.data = data

#     def dataframe(self) -> pd.DataFrame:
#         """Build the dataframe based on input format."""

#         return FACTORY_DICT[self.format].dataframe(self.data)
#         # return ExpandedDataFrame().dataframe(self.data)
