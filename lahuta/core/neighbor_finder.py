import numpy as np
from MDAnalysis.lib.nsgrid import FastNS

from lahuta.core.neighbors import NeighborPairs
from lahuta.utils.mda import mda_psuedobox_from_atomgroup


class NeighborSearch:
    """
    A class to handle atom related operations, such as preparing for computation, finding neighbours, etc.

    Attributes:
    ----------
    instance : Class Instance
        The instance of the class these methods originally belonged to.
    """

    def __init__(self, instance):
        """
        Initialize NeighborSearch class.

        Args:
        ----
        instance : Class Instance
            The instance of the class these methods originally belonged to.
        """
        self.instance = instance

    def compute(self, radius=5.0, ignore_hydrogens=True, res_dif=1):
        """
        Compute the neighbors of each atom in the Universe.

        Args:
        ----
        radius (float, optional): The cutoff radius. Default is 5.0.
        ignore_hydrogens (bool, optional): Whether to ignore hydrogens. Default is True.
        skip_adjacent (bool, optional): Whether to skip adjacent. Default is True.
        res_dif (int, optional): The residue difference to consider. Default is 1.

        Returns:
        -------
        NeighborPairs : np.ndarray
            An array of shape (n_atoms, n_neighbors) where each row contains the indices of the neighbors of the atom in the row.

        """
        if ignore_hydrogens:
            atomgroup = self.instance.universe.atoms.select_atoms("not name H*")
        else:
            atomgroup = self.instance.atoms

        pairs, distances = self.get_neighbors(atomgroup, radius)

        print("max pairs values: ", pairs.max(axis=0))

        if res_dif > 0:
            idx = self._remove_adjacent_residue_pairs(pairs, res_dif=res_dif)
            pairs = pairs[idx]
            distances = distances[idx]

        return NeighborPairs(self.instance, pairs, distances)

    def get_neighbors(self, atomgroup=None, radius=5.0):
        """
        Get the neighbors of an atomgroup.

        Args:
        ----
        atomgroup (AtomGroup): The atomgroup to get the neighbors of.

        Returns:
        ------
        tuple: A tuple containing two arrays. The first array of shape (n_pairs, 2) where each row contains the indices of the atoms in the pair. The second array of shape (n_pairs,) containing the distances of each pair.
        """
        if atomgroup is None:
            atomgroup = self.instance.atoms

        shift_coords, pseudobox = mda_psuedobox_from_atomgroup(atomgroup)

        gridsearch = FastNS(
            cutoff=radius, coords=shift_coords, box=pseudobox, pbc=False
        )
        neighbors = gridsearch.self_search()

        return atomgroup[neighbors.get_pairs()].indices, neighbors.get_pair_distances()

    def _remove_adjacent_residue_pairs(self, pairs, res_dif=1):
        """
        Remove pairs where the difference in residue ids is less than `res_dif`.

        Args:
        ----
        pairs (np.ndarray): An array of shape (n_pairs, 2) where each row is a pair of atom indices.
        res_dif (int, optional): The difference in residue ids to remove. Default is 1.

        Returns:
        -------
        np.ndarray: An array of shape (n_pairs, 2) where each row is a pair of atom indices.
        """
        resids = self.instance.universe.atoms.resids[pairs]
        return np.any(np.abs(resids - resids[:, ::-1]) > res_dif, axis=1)
