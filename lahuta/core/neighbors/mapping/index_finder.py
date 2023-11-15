"""Find indices of pairs based on a given query pairs."""
import logging
from collections import defaultdict, deque

import numpy as np
from numpy.typing import NDArray


class IndexFinder:
    """Helper class for finding indices of pairs based on a given query pairs.

    NumPy structured arrays are used to store the pairs of atom names, residue IDs, and residue names.
    The class provides the logic and API to find the indices of pairs that match a given query pairs.
    It supports emptry strings in the query pairs. Empty strings act as wildcards and match any value for the
    corresponding field in the pairs.

    It provides two methods to find the indices of pairs that match the query pairs:
    1. A fast lookup method that uses a dictionary to store the indices of pairs that match the query pairs. For this
         method to work, the query pairs must not contain empty strings in the 'resids' field.
    2. A slow method that uses a `deque` to store the indices of pairs that match the query pairs. This method works
            for any query pairs.

    Attributes:
        pairs (NDArray[np.void]): Array of pairs of atom names, residue IDs, and residue names.
        indices (list[int]): List of indices of pairs that match the query pairs.

    Methods:
        find_indices: Find the indices of pairs that match the query pairs.
    """

    def __init__(self, pairs: NDArray[np.void]):
        self.pairs = pairs
        self.indices: list[int] = []
        self._categorized_pairs = self._categorize()

    def _categorize(self) -> dict[tuple[str, str], list[int]]:
        categorized = defaultdict(list)
        for idx, (pair1, pair2) in enumerate(self.pairs):
            # Using 'resids' as a quick lookup key
            key = (pair1[1], pair2[1])
            categorized[key].append(idx)
        return categorized

    @staticmethod
    def _match_tuple(pair1: NDArray[np.void], pair2: NDArray[np.void]) -> bool:
        return all(x == y or x == "" or y == "" for x, y in zip(pair1, pair2, strict=True))

    @staticmethod
    def _has_empty_resids(pairs: NDArray[np.void]) -> np.bool_:
        return np.any(pairs["resids"] == "")

    def find_indices(self, ppairs: NDArray[np.void]) -> list[int]:
        """Find the indices of pairs that match the query pairs.

        Args:
            ppairs (NDArray[np.void]): Array of pairs of atom names, residue IDs, and residue names.

        Returns:
            list[int]: List of indices of pairs that match the query pairs.
        """
        if self._has_empty_resids(ppairs):
            logging.info("Found empty resids! Cannot use fast lookup. Using slow approach instead.")
            self._find_indices_using_deque(ppairs)
            return sorted(self.indices)

        self._find_indices_fast_lookup(ppairs)
        return self.indices

    def _find_indices_fast_lookup(self, ppairs: NDArray[np.void]) -> None:
        for pp1, pp2 in ppairs:
            resids_key = (pp1[1], pp2[1])
            if resids_key in self._categorized_pairs:
                candidates = self._categorized_pairs[resids_key]
                for i in candidates:
                    if self._match_tuple(self.pairs[i][0], pp1) and self._match_tuple(self.pairs[i][1], pp2):
                        self.indices.append(i)  # noqa: PERF401

    def _find_indices_using_deque(self, ppairs: NDArray[np.void]) -> None:
        candidates = deque(range(len(self.pairs)))

        for pp1, pp2 in ppairs:
            new_candidates: deque[int] = deque()
            for i in candidates:
                if self._match_tuple(self.pairs[i][0], pp1) and self._match_tuple(self.pairs[i][1], pp2):
                    self.indices.append(i)
                else:
                    new_candidates.append(i)
            candidates = new_candidates
