"""
Placeholder for the universe module.
"""

import lahuta.contacts as F
from lahuta.core.neighbors import NeighborPairs

from .base import ContactAnalysis


class VanDerWaalsContacts(ContactAnalysis):
    """A class to compute VdW contacts."""

    vdw_comp_factor = 0.1
    remove_clashes = True

    def compute(self) -> NeighborPairs:
        return F.vdw_neighbors(self.ns, self.vdw_comp_factor, self.remove_clashes)
