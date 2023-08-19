"""Handle atom related operations, including finding neighbors and preparation for computation.

This module contains the NeighborSearch class that handles atom related operations, 
including finding neighbors and preparation for computation.

Classes:
    NeighborSearch: Class to handle atom related operations, including finding neighbors 
                    and preparation for computation.

"""
from typing import Tuple

import numpy as np
from MDAnalysis.lib.nsgrid import FastNS
from numpy.typing import NDArray

from lahuta.lahuta_types.mdanalysis import AtomGroupType
from lahuta.utils.mda import mda_psuedobox_from_atomgroup

IndexPairs = NDArray[np.int32]
Distances = NDArray[np.float32]
PairsDistances = Tuple[IndexPairs, Distances]


class NeighborSearch:
    """Handle atom related operations, including finding neighbors and preparation for computation.

    The class provides methods to find neighbors of each atom in the universe and to remove pairs of atoms
    that are adjacent in the sequence.

    Args:
        mda (AtomGroupType): The AtomGroup containing the atoms.


    Attributes:
        ag_no_h (AtomGroup): Atom group of a universe excluding hydrogen atoms.
        og_resids (np.ndarray): The residue IDs of each atom in the universe.
    """

    def __init__(self, mda: AtomGroupType) -> None:
        mda_atoms, mda_universe = mda.atoms, mda.universe
        self.ag_no_h = mda_atoms.select_atoms("not name H*")
        self.og_resids = mda_universe.atoms.resids

    def compute(self, radius: float = 5.0, res_dif: int = 1) -> PairsDistances:
        """Compute the neighbors of each atom in the Universe.

        Args:
            radius (float, optional): The cutoff radius. Default is 5.0.
            res_dif (int, optional): The residue difference to consider. Default is 1.

        Returns:
            PairsDistances: A tuple containing the pairs of atom indices and the distances.
        """
        pairs, distances = self.get_neighbors(radius)

        if res_dif > 0:
            idx = self.remove_adjacent_residue_pairs(pairs, res_dif=res_dif)
            pairs = pairs[idx]
            distances = distances[idx]

        return pairs, distances

    def get_neighbors(self, radius: float) -> PairsDistances:
        """Get the neighbors of an AtomGroup.

        Args:
            radius (float, optional): The cutoff radius.

        Returns:
            PairsDistances: A tuple containing the pairs of atom indices and the distances.
        """
        # check for dimensions
        if not hasattr(self.ag_no_h.universe, "dimensions") or self.ag_no_h.universe.dimensions is None:
            positions, dimensions = mda_psuedobox_from_atomgroup(self.ag_no_h)
            pbc = False
        else:
            positions = self.ag_no_h.positions
            dimensions = self.ag_no_h.universe.dimensions
            pbc = True

        gridsearch = FastNS(cutoff=radius, coords=positions, box=dimensions, pbc=pbc)
        neighbors = gridsearch.self_search()

        return (
            self.ag_no_h[neighbors.get_pairs()].ix,
            neighbors.get_pair_distances(),
        )

    def remove_adjacent_residue_pairs(self, pairs: NDArray[np.int32], res_dif: int = 1) -> NDArray[np.bool_]:
        """Remove pairs where the difference in residue ids is less than `res_dif`.

        Args:
            pairs (NDArray[np.int32]): An array of shape (n_pairs, 2) where each row is a pair of atom indices.
            res_dif (int, optional): The difference in residue ids to remove. Default is 1.

        Returns:
            NDArray[np.bool_]: An array of shape (n_pairs,) containing the indices of the pairs to keep.
        """
        resids = self.og_resids[pairs]
        # return np.any(np.abs(resids - resids[:, ::-1]) > res_dif, axis=1) # noqa: ERA001
        mask: NDArray[np.bool_] = np.abs(np.diff(resids, axis=1)) > res_dif
        return np.ravel(mask)
