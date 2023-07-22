from lahuta.core.neighbors import NeighborPairs

from .atom_plane import DEFAULT_CONTACT_DISTS, carbon_pi
from .base import ContactAnalysis


class CarbonPi(ContactAnalysis):
    distance = DEFAULT_CONTACT_DISTS["carbon_pi"]
    cache = False

    def compute(self) -> NeighborPairs:
        """Compute aromatic contacts."""
        return carbon_pi(self.ns, cache=self.cache)
