"""Defines a class for computing aromatic contacts using a class-based approach.
The AromaticContacts class inherits from the base ContactAnalysis class and 
implements the `compute` method for aromatic contact computation.

Class:
    AromaticContacts(ContactAnalysis): Computes aromatic contacts.
                                       

Example:
    luni = Luni(...)
    ns = luni.compute_neighbors()

    ac = AromaticContacts(ns)
    print(ac.results)

"""

import lahuta.contacts as F
from lahuta.config.defaults import CONTACTS
from lahuta.core.neighbors import NeighborPairs

from .base import ContactAnalysis


class AromaticContacts(ContactAnalysis):
    """Handles the computation of aromatic contacts in a molecular structure.

    Aromatic contacts are computed based on the interactions between atoms in
    aromatic rings found in proteins and ligands. Aromatic interactions,
    commonly found in biological systems, are characterized by π-stacking
    (face-to-face), T-shaped or edge-to-face configurations, and cation-π
    interactions.

    !!! tip "Definition"
        1. Both atoms belong to an aromatic ring.
        2. The distance between these two atoms does not exceed a defined distance cutoff.

    Attributes:
        ns (NeighborPairs): A NeighborPairs object containing the atom neighbor relationships in the system.
        distance (float): The maximum distance to consider for contact. See `lahuta.config.defaults.CONTACTS` for
            default values.

    ??? example "Example"
        ``` py
        luni = Luni(...)
        ns = luni.compute_neighbors()

        ac = AromaticContacts(ns)
        print(ac.results)
        ```
    """

    distance = CONTACTS["aromatic"]["distance"]

    def compute(self) -> NeighborPairs:
        """Compute aromatic contacts based on the neighbor pairs.

        Returns:
            NeighborPairs: A NeighborPairs object containing only aromatic contacts.
        """
        return F.aromatic_neighbors(self.ns, self.distance)
