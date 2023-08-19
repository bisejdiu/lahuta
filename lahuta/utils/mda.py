"""Create a pseudo box for PDB files and the shifting of coordinates within that box."""

import numpy as np
from numpy.typing import NDArray

from lahuta.lahuta_types.mdanalysis import AtomGroupType


def mda_psuedobox_from_atomgroup(
    ag: AtomGroupType, cutoff: float = 5.0
) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """Construct a pseudo box and shifts atomic coordinates to fit within this box.

    This function takes an AtomGroup object and a cutoff distance as input, and it creates
    a pseudo box large enough to include all atoms plus the specified cutoff. It then shifts
    the coordinates of the atoms to fit within this box.

    Args:
        ag (AtomGroupType): An AtomGroup object representing the group of atoms for which \
                             the pseudo box is to be created.
        cutoff (float, optional): The minimum distance that should be allowed between the \
                                  edges of the pseudo box and the atoms. Default is 5.0.

    Returns:
        tuple[NDArray[np.float32], NDArray[np.float32]]: The shifted coordinates of the atoms \
                                                        and the dimensions of the pseudo box.
    """
    lmax: float = ag.positions.max(axis=0)
    lmin: float = ag.positions.min(axis=0)
    pseudobox = np.zeros(6, dtype=np.float32)
    lengths = 1.1 * (lmax - lmin)
    lengths = np.maximum(lengths, 2 * cutoff)
    pseudobox = np.zeros(6, dtype=np.float32)
    pseudobox[:3] = lengths
    pseudobox[3:] = 90.0
    shift_coords: NDArray[np.float32] = ag.positions.copy()
    shift_coords -= lmin
    return shift_coords, pseudobox
