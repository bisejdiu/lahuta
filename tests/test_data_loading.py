import warnings
from pathlib import Path

import numpy as np

from lahuta import Luni


def read_pdb(pdb_file: str) -> Luni:
    """Read a PDB file and return a Luni object."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Luni(pdb_file)


def test_read_pdb() -> None:
    path_obj = Path(__file__).parent / "data" / "1KX2.pdb"
    u = read_pdb(str(path_obj))

    assert u.arc is not None
    assert u.arc.atoms.ids.size == 1249
    assert u.arc.residues is not None
    assert np.unique(u.arc.residues.resids).size == 82
    ca = u.to("mda").select_atoms("name CA")
    assert ca.indices is not None
    assert ca.indices.size == 82 - 1
