import warnings
from pathlib import Path

import numpy as np

from lahuta.core.universe import Universe


def read_pdb(pdb_file: str) -> Universe:
    """Read a PDB file and return a Universe object."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Universe(pdb_file)


def test_read_pdb():

    path_obj = Path(__file__).parent / "data" / "1KX2.pdb"
    u = read_pdb(str(path_obj))

    assert u.atoms.ids.size == 1249
    assert u.residues is not None
    assert np.unique(u.residues.resids).size == 82
    ca = u.to("mda").select_atoms("name CA")
    assert ca.indices is not None
    assert ca.indices.size == 82 - 1
