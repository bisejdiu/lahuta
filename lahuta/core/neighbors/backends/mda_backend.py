"""Handle atom related operations, including finding neighbors and preparation for computation.

This module contains the MDAnalysisNeighborSearch class that handles atom related operations,
including finding neighbors and preparation for computation.

Classes:
    MDAnalysisNeighborSearch: Class to handle atom related operations, including finding neighbors
                    and preparation for computation.

"""

from typing import Literal, Optional

import numpy as np
from MDAnalysis.lib.nsgrid import FastNS
from numpy.typing import NDArray

from lahuta.core.neighbors import BaseNeighborSearch, PairsDistances
from lahuta.utils.mda import mda_psuedobox_from_atomgroup


class MDAnalysisNeighborSearch(BaseNeighborSearch):
    """Handle atom related operations, including finding neighbors and preparation for computation.

    The class provides methods to find neighbors of each atom in the universe and to remove pairs of atoms
    that are adjacent in the sequence.

    Args:
        mda (AtomGroupType): The AtomGroup containing the atoms.


    Attributes:
        ag_no_h (AtomGroup): Atom group of a universe excluding hydrogen atoms.
        og_resids (np.ndarray): The residue IDs of each atom in the universe.
    """

    def compute(
        self,
        radius: float = 5.0,
        res_dif: int = 1,
        chain_type: Optional[Literal["inter", "intra"]] = None,
        image: Optional[Literal["inter", "intra"]] = None,
    ) -> PairsDistances:
        """Compute the neighbors of each atom in the Universe.

        Args:
            radius (float, optional): The cutoff radius. Default is 5.0.
            res_dif (int, optional): The residue difference to consider. Default is 1.
            chain_type (Optional[Literal["inter", "intra"]], optional): The type of chain to keep. Default is None.
            image (Optional[Literal["inter", "intra"]], optional): The type of image to keep. Default is None.

        Returns:
            PairsDistances: A tuple containing the pairs of atom indices and the distances.
        """
        if image is not None:
            raise NotImplementedError("Image handling not implemented yet.")

        pairs, distances = self.get_neighbors(radius)
        print("(Lahuta DEBUG) mda FastNS: ", pairs.shape, distances.shape)

        if res_dif > 0:
            idx = self.remove_adjacent_residue_pairs(pairs, res_dif=res_dif)
            pairs = pairs[idx]
            distances = distances[idx]

        if chain_type is not None:
            idx = self.filter_chain_type(pairs, chain_type=chain_type)
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
            positions, dimensions = mda_psuedobox_from_atomgroup(self.ag_no_h, cutoff=radius)
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
        # return np.any(np.abs(resids - resids[:, ::-1]) > res_dif, axis=1)
        mask: NDArray[np.bool_] = np.abs(np.diff(resids, axis=1)) > res_dif
        return np.ravel(mask)

    def filter_chain_type(self, pairs: NDArray[np.int32], chain_type: Literal["inter", "intra"]) -> NDArray[np.bool_]:
        """Remove pairs where the chain ids are the same.

        Args:
            pairs (NDArray[np.int32]): An array of shape (n_pairs, 2) where each row is a pair of atom indices.
            chain_type (Literal["inter", "intra"]): The type of chain to keep.

        Returns:
            NDArray[np.bool_]: An array of shape (n_pairs,) containing the indices of the pairs to keep.
        """
        assert chain_type in ["inter", "intra"], f"Invalid chain_type: {chain_type}. Must be 'inter' or 'intra'."

        chain_ids: NDArray[np.int32] = self.chain_ids[pairs]
        if chain_type == "inter":
            return chain_ids[:, 0] != chain_ids[:, 1]  # type: ignore

        return chain_ids[:, 0] == chain_ids[:, 1]  # type: ignore
