import warnings
from pathlib import Path

import numpy as np

from lahuta.core.universe import Universe


def read_pdb(pdb_file: Path) -> Universe:
    """Read a PDB file and return a Universe object."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Universe(pdb_file)


def test_read_pdb():

    u = read_pdb(Path(__file__).parent / "data" / "1KX2.pdb")

    assert u.atoms.indices.size == 1249
    assert u.residues is not None
    assert len(u.residues.indices) == 82
    ca = u.select_atoms("name CA")
    assert ca.indices is not None
    assert ca.indices.size == 82 - 1
