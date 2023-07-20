"""
Placeholder for the universe module.
"""

import lahuta.contacts as F

from .base import ContactAnalysis


class VanDerWaalsContacts(ContactAnalysis):
    """A class to compute VdW contacts."""

    vdw_comp_factor = 0.1
    remove_clashes = True

    def compute(self):
        return F.vdw_neighbors(self.ns, self.vdw_comp_factor, self.remove_clashes)
