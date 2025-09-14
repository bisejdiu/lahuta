"""Utility functions for MSA."""
import numpy as np
from numpy.typing import NDArray


def direct_merge(
    prot_resindices: NDArray[np.int32],
    nonprot_resindices: NDArray[np.int32],
    mapped_prot_resindices: NDArray[np.int32],
    mapped_nonprot_resindices: NDArray[np.int32],
) -> NDArray[np.int32]:
    """Sort mapped residue indices based on the original residue indices.

    Args:
        prot_resindices (NDArray[np.int32]): Array of protein residue indices.
        nonprot_resindices (NDArray[np.int32]): Array of non-protein residue indices.
        mapped_prot_resindices (NDArray[np.int32]): Array of mapped protein residue indices.
        mapped_nonprot_resindices (NDArray[np.int32]): Array of mapped non-protein residue indices.

    Returns:
        NDArray[np.int32]: Sorted array of mapped residue indices.
    """
    n = len(prot_resindices)
    m = len(nonprot_resindices)
    mapped_resindices = np.empty(n + m, dtype=mapped_prot_resindices.dtype)

    i = j = k = 0
    while i < n and j < m:
        if prot_resindices[i] < nonprot_resindices[j]:
            mapped_resindices[k] = mapped_prot_resindices[i]
            i += 1
        else:
            mapped_resindices[k] = mapped_nonprot_resindices[j]
            j += 1
        k += 1

    while i < n:
        mapped_resindices[k] = mapped_prot_resindices[i]
        i += 1
        k += 1

    while j < m:
        mapped_resindices[k] = mapped_nonprot_resindices[j]
        j += 1
        k += 1

    return mapped_resindices
