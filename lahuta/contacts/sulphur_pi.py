from .atom_plane import DEFAULT_CONTACT_DISTS, sulphur_pi
from .base import ContactAnalysis


class SulphurPi(ContactAnalysis):
    distance = DEFAULT_CONTACT_DISTS["sulphur_pi"]
    cache = False

    def compute(self):
        """Compute aromatic contacts."""
        return sulphur_pi(self.ns)
