import warnings

import numpy as np

from lahuta import Luni
from lahuta.tests import X2


def read_pdb(pdb_file: str) -> Luni:
    """Read a PDB file and return a Luni object."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Luni(pdb_file)


def test_read_pdb() -> None:
    u = read_pdb(str(X2()))

    assert u.arc is not None
    assert u.arc.atoms.ids.size == 1249
    assert u.arc.residues is not None
    assert np.unique(u.arc.residues.resids).size == 82
    ca = u.to("mda").select_atoms("name CA")
    assert ca.indices is not None
    assert ca.indices.size == 82 - 1
