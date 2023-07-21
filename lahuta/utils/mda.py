from typing import Tuple

import numpy as np
from numpy.typing import NDArray

from lahuta.types.mdanalysis import AtomGroupType


def mda_psuedobox_from_atomgroup(
    ag: AtomGroupType, cutoff: float = 5.0
) -> Tuple[NDArray[np.float_], NDArray[np.float_]]:
    """
    Create a psuedobox for MDAnalysis to use with PDB files.
    Also shifts the coordinates so they are within the box definition.
    """
    lmax = ag.positions.max(axis=0)
    lmin = ag.positions.min(axis=0)
    pseudobox = np.zeros(6, dtype=np.float32)
    lengths = 1.1 * (lmax - lmin)
    lengths = np.maximum(lengths, 2 * cutoff)
    pseudobox = np.zeros(6, dtype=np.float32)
    pseudobox[:3] = lengths
    pseudobox[3:] = 90.0
    shift_coords = ag.positions.copy()
    shift_coords -= lmin
    return shift_coords, pseudobox
