import numpy as np


def mda_psuedobox_from_atomgroup(atomgroup):
    """
    Create a psuedobox for MDAnalysis to use with PDB files.
    Also shifts the coordinates so they are within the box definition.
    """
    lmax = atomgroup.positions.max(axis=0)
    lmin = atomgroup.positions.min(axis=0)
    pseudobox = np.zeros(6, dtype=np.float32)
    pseudobox[:3] = 5 * (lmax - lmin)
    pseudobox[3:] = 90.0
    shift_coords = atomgroup.positions.copy()
    shift_coords -= lmin
    return shift_coords, pseudobox
