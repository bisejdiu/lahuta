"""
Placeholder for the universe module.
"""

from typing import Union

import numpy as np
from lahuta.config.defaults import CONTACTS
from lahuta.utils.array_utils import intersection
from MDAnalysis.lib import distances as mda_distances
from openbabel import openbabel as ob

from ..core.groups import AtomGroup
from ..core.neighbors import NeighborPairs
from ..core.universe import Universe
from ..utils.array_utils import non_matching_indices
from .protocol import ContactBase


class AtomPlaneContacts:
    def __init__(self, ua: Union[Universe, AtomGroup]):
        self.ua = ua

    def compute_contacts(self, **kwargs):

        rings = perceive_rings(self.ua.universe.mol)
        reference = np.array([ring["center"] for ring in rings])
        max_cutoff = CONTACTS["aromatic"]["met_sulphur_aromatic_distance"]
        pairs, d = get_mda_neighbors(self.ua, reference, max_cutoff=max_cutoff)
        processed_pairs = self.ua.atoms[pairs].indices

        atomgroup = self.ua.universe.select_atoms("not element H")
        processed_pairs = atomgroup[pairs].indices

        neighbors = NeighborPairs(self.ua.atoms, processed_pairs, d)
        no_aroms = neighbors.difference(neighbors.type_filter("aromatic", col=1))

        ix = intersection(processed_pairs, no_aroms.pairs)
        pix = pairs[ix]

        ring_centers = np.array([ring["center"] for ring in rings])
        ring_normals = np.array([ring["normal"] for ring in rings])

        ring_centers = np.take(ring_centers, pix[:, 0], axis=0)
        ring_normals = np.take(ring_normals, pix[:, 0], axis=0)

        angles = vector_angle(
            ring_normals, ring_centers - self.ua.atoms[pix[:, 1]].positions
        )

        result = []
        for (ix1, ix2), angle in zip(processed_pairs, angles):
            result.append([ix1, ix2, angle])

        return sorted(result, key=lambda x: x[2])


def get_mda_neighbors(ua, reference, max_cutoff=5.0):
    ag = ua.select_atoms("not element H")

    pairs, distances = mda_distances.capped_distance(
        reference, ag.positions, max_cutoff, return_distances=True
    )

    return pairs, distances


def vector_angle(v1, v2):

    dot = np.einsum("ij,ij->i", v1, v2 / np.linalg.norm(v2, axis=1)[..., np.newaxis])
    angle = np.sign(dot) * np.arccos(dot)

    angle = np.where(angle < 0, angle + np.pi, angle)

    return np.degrees(angle)


def perceive_rings(mol):
    o = []

    for e, ob_ring in enumerate(mol.GetSSSR()):

        if not ob_ring.IsAromatic():
            continue

        center = ob.vector3()
        normal = ob.vector3()

        ob_ring.findCenterAndNormal(center, normal, ob.vector3())

        # CONVERT CENTER AND NORMALS TO NUMPY
        center = np.array([center.GetX(), center.GetY(), center.GetZ()])
        normal = np.array([normal.GetX(), normal.GetY(), normal.GetZ()])

        o.append(
            {
                "ring_id": e,
                "center": center,
                "normal": normal,
                "atoms": sorted([atom for atom in ob_ring._path]),
            }
        )

    return o
