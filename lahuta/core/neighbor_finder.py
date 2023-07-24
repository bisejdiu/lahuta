from typing import Tuple

import numpy as np
from MDAnalysis.lib.nsgrid import FastNS  # type: ignore
from numpy.typing import NDArray

from lahuta.lahuta_types.mdanalysis import AtomGroupType
from lahuta.utils.mda import mda_psuedobox_from_atomgroup


class NeighborSearch:
    """
    A class to handle atom related operations, such as preparing for computation, finding neighbours, etc.

    Attributes:
    ----------
    instance : Class Instance
        The instance of the class these methods originally belonged to.
    """

    def __init__(self, mda: AtomGroupType) -> None:
        """
        Initialize NeighborSearch class.

        Args:
        ----
        instance : Class Instance
            The instance of the class these methods originally belonged to.
        """
        # mda = instance.to("mda")  # .copy()
        mda_atoms, mda_universe = mda.atoms, mda.universe
        self.ag_no_h = mda_atoms.select_atoms("not name H*")
        self.og_resids = mda_universe.atoms.resids

        # self.instance = instance

    def compute(self, radius: float = 5.0, res_dif: int = 1) -> Tuple[NDArray[np.int32], NDArray[np.float32]]:
        """
        Compute the neighbors of each atom in the Universe.

        Args:
        ----
        radius (float, optional): The cutoff radius. Default is 5.0.
        skip_adjacent (bool, optional): Whether to skip adjacent. Default is True.
        res_dif (int, optional): The residue difference to consider. Default is 1.

        Returns:
        -------
        NeighborPairs : np.ndarray
            An array of shape (n_atoms, n_neighbors) where each row contains the indices of the neighbors of the atom in the row.

        """

        pairs, distances = self.get_neighbors(radius)

        if res_dif > 0:
            idx = self._remove_adjacent_residue_pairs(pairs, res_dif=res_dif)
            pairs = pairs[idx]
            distances = distances[idx]
        # else:
        #     pairs = og_pairs

        # assert pairs.shape[0] != 0, "No neighbors found. Try increasing the radius."
        return pairs, distances

        # return NeighborPairs(self.instance, pairs, distances)

    def get_neighbors(self, radius: float = 5.0) -> Tuple[NDArray[np.int32], NDArray[np.float32]]:
        """
        Get the neighbors of an atomgroup.

        Args:
        ----
        atomgroup (AtomGroup): The atomgroup to get the neighbors of.

        Returns:
        ------
        tuple: A tuple containing two arrays. The first array of shape (n_pairs, 2) where each row contains the indices of the atoms in the pair. The second array of shape (n_pairs,) containing the distances of each pair.
        """
        # def fast_ns_wrapper(atomgroup, radius):
        #     shift_coords, pseudobox = mda_psuedobox_from_atomgroup(atomgroup)
        #     gridsearch = FastNS(
        #         cutoff=radius, coords=shift_coords, box=pseudobox, pbc=False
        #     )
        #     neighbors = gridsearch.self_search()
        #     return neighbors

        shift_coords, pseudobox = mda_psuedobox_from_atomgroup(self.ag_no_h)

        # TODO: handle pbc
        gridsearch = FastNS(cutoff=radius, coords=shift_coords, box=pseudobox, pbc=False)  # type: ignore
        neighbors = gridsearch.self_search()  # type: ignore

        return (
            self.ag_no_h[neighbors.get_pairs()].indices,  # type: ignore
            neighbors.get_pair_distances(),  # type: ignore
        )

    def _remove_adjacent_residue_pairs(self, pairs: NDArray[np.int32], res_dif: int = 1) -> NDArray[np.bool_]:
        """
        Remove pairs where the difference in residue ids is less than `res_dif`.

        Args:
        ----
        pairs (np.ndarray): An array of shape (n_pairs, 2) where each row is a pair of atom indices.
        res_dif (int, optional): The difference in residue ids to remove. Default is 1.

        Returns:
        -------
        np.ndarray: An array of shape (n_pairs,) containing the indices of the pairs to keep.
        """

        resids = self.og_resids[pairs]
        # return np.any(np.abs(resids - resids[:, ::-1]) > res_dif, axis=1)
        mask = np.abs(np.diff(resids, axis=1)) > res_dif
        return np.ravel(mask)
