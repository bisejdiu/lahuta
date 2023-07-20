from .atom_plane import DEFAULT_CONTACT_DISTS, donor_pi
from .base import ContactAnalysis


class DonorPi(ContactAnalysis):
    distance = DEFAULT_CONTACT_DISTS["donor_pi"]
    cache = False

    def compute(self):
        """Compute aromatic contacts."""
        return donor_pi(self.ns)
