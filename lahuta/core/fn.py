"""Handle atom related operations, including finding neighbors and preparation for computation."""
from concurrent.futures import ThreadPoolExecutor
from typing import Literal, Optional

import gemmi
import numpy as np
from gemmi import ContactSearch, NeighborSearch
from numpy.typing import NDArray

from lahuta.lahuta_types.gemmi import SearchResults, Structure
from lahuta.lahuta_types.mdanalysis import AtomGroupType


class GemmiNeighbors:
    """Handle atom related operations, including finding neighbors and preparation for computation.

    The class provides methods to find neighbors of each atom in the universe and to remove pairs of atoms
    that are adjacent in the sequence.

    Args:
        mda (AtomGroupType): The AtomGroup containing the atoms.
        structure (Structure): The Structure containing the atoms.
        radius (float, optional): The cutoff radius. Default is 5.0.

    Attributes:
        ag_no_h (AtomGroup): Atom group of a universe excluding hydrogen atoms.
        og_resids (np.ndarray): The residue IDs of each atom in the universe.
    """

    def __init__(self, mda: AtomGroupType, structure: Structure, radius: float=5.0):
        self.ag_no_h = mda.select_atoms("not name H*")
        self.og_resids = mda.universe.atoms.resids
        self.chain_ids = mda.universe.atoms.chainIDs

        structure.assign_serial_numbers()
        self.structure = structure
        self.radius = radius
        self.results: list[SearchResults] = []
        
    def _get_contacts(self, radius: float) -> list[SearchResults]:
        ns = NeighborSearch(self.structure[0], self.structure.cell, radius)
        ns.populate(include_h=False)
        cs = ContactSearch(radius)
        cs.ignore = gemmi.ContactSearch.Ignore.Nothing

        return cs.find_contacts(ns) # type: ignore

    def compute(
            self,
            radius: float = 5.0,
            res_dif: int = 1, 
            n_threads: int = 1, 
            chain_type: Optional[Literal["inter", "intra"]] = None, 
            image: Optional[Literal["inter", "intra"]] = None
    ) -> tuple[NDArray[np.int32], NDArray[np.float32]]:
        """Compute the neighbors of each atom in the loaded system.

        Args:
            radius (float, optional): The cutoff radius. Default is 5.0.
            res_dif (int, optional): The residue difference to consider. Default is 1.
            n_threads (int, optional): The number of threads to use. Default is 1.
            chain_type (Literal["inter", "intra"], optional): The type of chain to keep. Default is None (keep all).
            image (Literal["inter", "intra"], optional): The type of image to keep. Default is None (keep all).

        Returns:
            tuple[NDArray[np.int32], NDArray[np.float32]]: A tuple containing the pairs of atom indices and the distances.
        """
        results = self._get_contacts(radius)
        pairs, distances, image_ids = self._compute(results, n_threads)

        if res_dif > 0:
            idx = self.remove_adjacent_residue_pairs(pairs, res_dif=res_dif)
            pairs = pairs[idx]
            distances = distances[idx]

        if chain_type is not None:
            idx = self.filter_chain_type(pairs, chain_type=chain_type)
            pairs = pairs[idx]
            distances = distances[idx]

        if image is not None:
            idx = self.filter_image(image_ids, pairs, image=image)
            pairs = pairs[idx]
            distances = distances[idx]

        return pairs, distances

    def _compute(
            self, 
            results: list[SearchResults], 
            n_threads: int
            ) -> tuple[NDArray[np.int32], NDArray[np.float32], NDArray[np.int32]]:
        n_results = len(results)
        pairs = np.empty((n_results, 2), dtype=np.int32)
        distances = np.empty(n_results, dtype=np.float32)
        image_ids = np.empty(n_results, dtype=np.int32)

        slice_size = n_results // n_threads

        with ThreadPoolExecutor() as executor:
            futures = {executor.submit(self.fetch_attributes, results[i:i+slice_size]): i for i in range(0, n_results, slice_size)}

        # Collate results
        for future, start_idx in futures.items():
            slice_pairs, slice_distances, slice_image_ids = future.result()
            end_idx = start_idx + len(slice_distances)
            pairs[start_idx:end_idx, :] = slice_pairs
            distances[start_idx:end_idx] = slice_distances
            image_ids[start_idx:end_idx] = slice_image_ids

        pairs -= 1 

        return pairs, distances, image_ids

    @staticmethod
    def fetch_attributes(
        results_slice: list[SearchResults]
        ) -> tuple[NDArray[np.int32], NDArray[np.float32], NDArray[np.int32]]:
        """Fetch the attributes from a slice of SearchResults.

        Args:
            results_slice (list[SearchResults]): A slice of SearchResults.

        Returns:
            tuple[NDArray[np.int32], NDArray[np.float32], NDArray[np.int32]]: A tuple containing the pairs of atom indices, the distances, and the image ids.
        """
        n = len(results_slice)
        slice_pairs = np.empty((n, 2), dtype=np.int32)
        slice_distances = np.empty(n, dtype=np.float32)
        slice_image_ids = np.empty(n, dtype=np.int32)

        for i, r in enumerate(results_slice):
            slice_pairs[i, 0] = r.partner1.atom.serial
            slice_pairs[i, 1] = r.partner2.atom.serial
            slice_distances[i] = r.dist
            slice_image_ids[i] = r.image_idx

        return slice_pairs, slice_distances, slice_image_ids

    def remove_adjacent_residue_pairs(self, pairs: NDArray[np.int32], res_dif: int = 1) -> NDArray[np.bool_]:
        """Remove pairs where the difference in residue ids is less than `res_dif`.

        Args:
            pairs (NDArray[np.int32]): An array of shape (n_pairs, 2) where each row is a pair of atom indices.
            res_dif (int, optional): The difference in residue ids to remove. Default is 1.

        Returns:
            NDArray[np.bool_]: An array of shape (n_pairs,) containing the indices of the pairs to keep.
        """
        resids = self.og_resids[pairs]
        mask = np.abs(np.diff(resids, axis=1)) > res_dif
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
            return chain_ids[:, 0] != chain_ids[:, 1] # type: ignore
        
        return chain_ids[:, 0] == chain_ids[:, 1] # type: ignore
    
    def filter_image(
            self, 
            image_ids: NDArray[np.int32], 
            pairs: NDArray[np.int32], 
            image: Literal["inter", "intra"]
        ) -> NDArray[np.bool_]:
        """Remove pairs where the image ids are the same.

        Args:
            image_ids (NDArray[np.int32]): An array of shape (n_pairs,) containing the image ids.
            pairs (NDArray[np.int32]): An array of shape (n_pairs, 2) where each row is a pair of atom indices.
            image (Literal["inter", "intra"]): The type of chain to keep.

        Returns:
            NDArray[np.bool_]: An array of shape (n_pairs,) containing the indices of the pairs to keep.
        """
        assert image in ["inter", "intra"], f"Invalid image: {image}. Must be 'inter' or 'intra'."
        image_pairs = image_ids[pairs]
        if image == "inter":
            return image_pairs[:, 0] != image_pairs[:, 1] # type: ignore

        return image_pairs[:, 0] == image_pairs[:, 1] # type: ignore
