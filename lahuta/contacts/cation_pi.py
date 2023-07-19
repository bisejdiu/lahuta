from .atom_plane import DEFAULT_CONTACT_DISTS, cation_pi
from .base import ContactAnalysis


class CationPi(ContactAnalysis):
    distance = DEFAULT_CONTACT_DISTS["cation_pi"]
    cache = False

    def compute(self):
        """Compute aromatic contacts."""
        return cation_pi(self.ns)
