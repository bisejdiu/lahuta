"""
Placeholder for the universe module.
"""

import itertools
from dataclasses import dataclass
from typing import Any, Dict, List, Literal, Optional, Protocol, Union

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib import distances as mda_distances
from openbabel import openbabel as ob

from lahuta.config.defaults import CONTACTS
from lahuta.contacts.protocol import ContactBase
from lahuta.core.universe import Universe

from ..core.neighbors import NeighborPairs
from ..utils.writers import FACTORY_DICT


class PlanePlaneContacts:
    def __init__(self, ua: Union[Universe, mda.AtomGroup]):
        self.ua = ua
        self.rings = perceive_rings(self.ua.to("mol"))

    def compute_contacts(self, **kwargs):
        ring_ids = np.arange(len(self.rings))

        reference = np.array([x["center"] for x in self.rings])
        distances = mda_distances.self_distance_array(reference, box=None)
        differences = reference[:, None] - reference[None, :]

        # results = []
        pairs = []
        pair_distances = []
        pair_contact_label = []
        ring_atom_indices = []
        theta_angles = []
        normal_angles = []

        ring_pairs = itertools.combinations(ring_ids, 2)
        for ix, (ix1, ix2) in enumerate(ring_pairs):
            ring1 = self.rings[ix1]
            ring2 = self.rings[ix2]

            if distances[ix] > CONTACTS["aromatic"]["centroid_distance"]:
                continue

            theta_point = differences[ix1, ix2]
            normal_angle = vector_angle_1d(ring1["normal"], ring2["normal"])
            theta = vector_angle_1d(ring1["normal"], theta_point)

            int_type = assign_pp_contact_type(normal_angle, theta)

            if int_type is None:
                continue

            resid1 = self.ua.to("mda").atoms[ring1["atoms"]].indices[0]
            resid2 = self.ua.to("mda").atoms[ring2["atoms"]].indices[0]

            pairs.append([resid1, resid2])
            pair_distances.append(distances[ix])
            pair_contact_label.append(int_type)
            ring_atom_indices.append([ring1["atoms"], ring2["atoms"]])
            theta_angles.append(theta)
            normal_angles.append(normal_angle)

            # results.append([ring1, ring2, distances[ix], normal_angle, theta, int_type])

        self.pairs = np.array(pairs)
        self.distances = np.array(pair_distances)
        self.contact_labels = np.array(pair_contact_label)
        self.ring_atom_indices = ring_atom_indices
        self.theta_angles = np.array(theta_angles)
        self.normal_angles = np.array(normal_angles)
        # return results


def assign_pp_contact_type(normal_angle, theta):
    """Assigns a contact type based on the normal angle and theta angle. Avoid using if/else statements."""

    # define the different types
    types = ["FF", "OF", "EE", "FT", "OT", "ET", "FE", "OE", "EF"]

    # define the different ranges
    ranges = [
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

    # define the different conditions
    conditions = [
        (normal_angle >= r[0])
        and (normal_angle <= r[1])
        and (theta >= r[2])
        and (theta <= r[3])
        for r in ranges
    ]

    # get the index of the first condition that is True
    index = next(i for i, x in enumerate(conditions) if x)

    # return the type
    return types[index]


def vector_angle(v1, v2):
    """Calculate the angle between two vectors.

    Parameters
    ----------
    v1 : np.ndarray
        First vector. It must be already normalized.
    v2 : np.ndarray
        Second vector.

    Returns
    -------
    float
        Angle between vectors in degrees.
    """

    dot = np.einsum("ij,ij->i", v1, v2 / np.linalg.norm(v2, axis=1)[..., np.newaxis])
    angle = np.sign(dot) * np.arccos(dot)

    angle = np.where(angle < 0, angle + np.pi, angle)

    return np.degrees(angle)


def vector_angle_1d(v1, v2):
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    dot = np.dot(v1_u, v2_u)
    angle = np.arccos(np.clip(dot, -1.0, 1.0))
    angle = np.sign(dot) * angle
    angle = np.degrees(angle)

    if angle < 0:
        angle = 180 + angle

    return angle


def perceive_rings(mol):
    aromatics = []

    for _, ob_ring in enumerate(mol.GetSSSR()):
        if not ob_ring.IsAromatic():
            continue

        center = ob.vector3()
        normal = ob.vector3()

        ob_ring.findCenterAndNormal(center, normal, ob.vector3())

        center = np.array([center.GetX(), center.GetY(), center.GetZ()])
        normal = np.array([normal.GetX(), normal.GetY(), normal.GetZ()])
        atoms = sorted([atom for atom in ob_ring._path])
        aromatics.append(
            {
                "atoms": atoms,
                "center": center,
                "normal": normal,
            }
        )

    return aromatics


def process_ring_information(indices, rings):
    """Process ring information for a given set of indices."""
    ring_atoms = {ix: ring["atoms"] for ix, ring in enumerate(rings) if ix in indices}
    ring_atoms = [ring_atoms[ix] for ix in indices]
    first_ring_atoms = [ring[0] for ring in ring_atoms]

    return first_ring_atoms, ring_atoms


class APDataFrameFactory:
    """A class for storing and manipulating contact data."""

    def __init__(
        self,
        apcontacts: Any,
        rings: Any,
        label: str = "ring",
        df_format: Literal["print", "compact", "expanded"] = "print",
    ):
        """A class for storing and manipulating contact data."""

        nag = apcontacts._atoms[apcontacts.pairs]

        first_rings, ring_atoms = process_ring_information(nag[:, 0].indices, rings)

        col1, col2 = nag.universe.atoms[first_rings], nag[:, 1]

        self.methods = ["resids", "resnames", "names", "indices"]
        self.format = df_format

        data = {}
        for method in self.methods:
            data[f"residue1_{method}"] = getattr(col1, method)
            data[f"residue2_{method}"] = getattr(col2, method)

        data["ring_atoms"] = ring_atoms
        data["distances"] = apcontacts.distances
        data["labels"] = [label] * len(apcontacts.distances)

        self.data = data

    def dataframe(self) -> pd.DataFrame:
        """Build the dataframe based on input format."""

        return FACTORY_DICT[self.format].dataframe(self.data)
        # return ExpandedDataFrame().dataframe(self.data)


class PPDataFrameFactory:
    """A class for storing and manipulating contact data."""

    def __init__(
        self,
        pcontacts: Any,
        df_format: Literal["print", "compact", "expanded"] = "print",
    ):
        """A class for storing and manipulating contact data."""

        nag = pcontacts.ua.to("mda").atoms[pcontacts.pairs]
        col1, col2 = nag[:, 0], nag[:, 1]

        self.methods = ["resids", "resnames", "names", "indices"]
        self.format = df_format

        data = {}
        for method in self.methods:
            data[f"residue1_{method}"] = getattr(col1, method)
            data[f"residue2_{method}"] = getattr(col2, method)

        data["ring1_atoms"] = [ring[0] for ring in pcontacts.ring_atom_indices]
        data["ring2_atoms"] = [ring[1] for ring in pcontacts.ring_atom_indices]
        data["distances"] = pcontacts.distances
        data["theta_angles"] = pcontacts.theta_angles
        data["normal_angles"] = pcontacts.normal_angles
        data["labels"] = pcontacts.contact_labels

        self.data = data

    def dataframe(self) -> pd.DataFrame:
        """Build the dataframe based on input format."""

        return FACTORY_DICT[self.format].dataframe(self.data)
        # return ExpandedDataFrame().dataframe(self.data)
