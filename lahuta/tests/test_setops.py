import warnings
from pathlib import Path

import numpy as np

from lahuta.core.neighbors import NeighborPairs
from lahuta.core.universe import Universe


def test_intersection(repeat: int = 10):
    """Test the intersection of two NeighborPairs objects."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = Universe(Path(__file__).parent / "data" / "1KX2.pdb")
        n = u.compute_neighbors()

    for _ in range(repeat):
        sample_size = np.random.randint(1, n.pairs.shape[0])
        i = np.random.choice(n.pairs.shape[0], sample_size)
        i = np.unique(np.sort(i))
        m = NeighborPairs(u.atoms, n.pairs[i], n.distances[i])

        assert (
            i.size
            == m.pairs.shape[0]
            == n.pairs[i].shape[0]
            == n.distances[i].shape[0]
            == m.distances.shape[0]
        )
