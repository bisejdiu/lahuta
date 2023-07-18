"""
Placeholder for the universe module.
"""

import lahuta.contacts as F
from lahuta.config.defaults import CONTACTS

# from ..config.atoms import METALS
from .base import ContactAnalysis


class MetalicContacts(ContactAnalysis):
    """A class to compute metal contacts."""

    distance = CONTACTS["metal"]["distance"]

    # def __init__(self, ns):
    #     self.metal_indices = (
    #         ns.atoms[ns.indices].select_atoms("element " + " ".join(METALS)).indices
    #     )

    #     super().__init__(ns)

    def compute(self):
        """Compute metal contacts."""
        return F.metalic_neighbors(self.ns, self.distance)
