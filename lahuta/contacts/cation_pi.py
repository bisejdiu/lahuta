from lahuta.core.neighbors import NeighborPairs

from .atom_plane import DEFAULT_CONTACT_DISTS, cation_pi
from .base import ContactAnalysis


class CationPi(ContactAnalysis):
    distance = DEFAULT_CONTACT_DISTS["cation_pi"]
    cache = False

    def compute(self) -> NeighborPairs:
        """Compute aromatic contacts."""
        return cation_pi(self.ns, cache=self.cache)
