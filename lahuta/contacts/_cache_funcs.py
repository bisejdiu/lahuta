from typing import Callable, Tuple

import numpy as np
from joblib import Memory
from MDAnalysis.lib import distances as mda_distances
from numpy.typing import NDArray

from lahuta.config.defaults import CONTACTS
from lahuta.core.neighbors import NeighborPairs
from lahuta.lahuta_types.mda_commands import CappedDistance, DistanceType
from lahuta.lahuta_types.mdanalysis import AtomGroupType
from lahuta.utils.math import calc_vec_line_angles


def compute_neighbors(
    positions: NDArray[np.float32], reference: NDArray[np.float32]
) -> Tuple[NDArray[np.int32], NDArray[np.float32]]:
    """Compute the neighbors between the reference and positions."""
    max_cutoff = CONTACTS["aromatic"]["met_sulphur_aromatic_distance"]

    wrapper: DistanceType = CappedDistance(mda_distances)
    pairs, distances = wrapper.capped_distance(reference, positions, max_cutoff, return_distances=True)
    return pairs, distances


def calc_ringnormal_pos_angle(
    ns: NeighborPairs, uv_atoms: AtomGroupType, ring_centers: NDArray[np.float32], ring_normals: NDArray[np.float32]
) -> NDArray[np.float32]:
    """Calculate the angle between the ring normal and the vector connecting the ring center and the atom."""
    selected_ring_centers = ring_centers[ns.pairs[:, 0]]
    selected_ring_normals = ring_normals[ns.pairs[:, 0]]

    return calc_vec_line_angles(
        selected_ring_normals,
        selected_ring_centers - uv_atoms[ns.pairs[:, 1]].positions,
    )


memory = Memory("cachedir", verbose=0)
compute_neighbors_cached: Callable[
    [NDArray[np.float32], NDArray[np.float32]], Tuple[NDArray[np.int32], NDArray[np.float32]]
]
compute_neighbors_cached = memory.cache(compute_neighbors)  # type: ignore
calc_ringnormal_pos_angle_cached: Callable[
    [NeighborPairs, AtomGroupType, NDArray[np.float32], NDArray[np.float32]], NDArray[np.float32]
]
calc_ringnormal_pos_angle_cached = memory.cache(calc_ringnormal_pos_angle)  # type: ignore
