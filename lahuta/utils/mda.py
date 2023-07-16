import numpy as np


def mda_psuedobox_from_atomgroup(atomgroup, cutoff=5.0):
    """
    Create a psuedobox for MDAnalysis to use with PDB files.
    Also shifts the coordinates so they are within the box definition.
    """
    lmax = atomgroup.positions.max(axis=0)
    lmin = atomgroup.positions.min(axis=0)
    pseudobox = np.zeros(6, dtype=np.float32)
    lengths = 1.1 * (lmax - lmin)
    lengths = np.maximum(lengths, 2 * cutoff)
    pseudobox = np.zeros(6, dtype=np.float32)
    pseudobox[:3] = lengths
    pseudobox[3:] = 90.0
    shift_coords = atomgroup.positions.copy()
    shift_coords -= lmin
    return shift_coords, pseudobox
