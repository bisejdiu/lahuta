"""Defines a class for computing ionic contacts using a class-based approach.
The IonicContacts class inherits from the base ContactAnalysis class and 
implements the `compute` method for ionic contact computation.

Class:
    IonicContacts(ContactAnalysis): Computes ionic contacts.
                                       

Example:
    luni = Luni(...)
    ns = luni.compute_neighbors()

    ionic = IonicContacts(ns)
    print(ionic.results)

"""

import lahuta.contacts as F
from lahuta.config.defaults import CONTACTS
from lahuta.core.neighbors import NeighborPairs

from .base import ContactAnalysis


class IonicContacts(ContactAnalysis):
    """Handle the computation of ionic contacts in a molecular system.

    Ionic contacts refer to the interactions between positively and negatively ionizable atoms,
    forming one of the primary types of electrostatic interactions.

    !!! tip "Definition"
        1. One atom must be positively ionizable.
        2. The other atom must be negatively ionizable.
        3. The distance between these two atoms does not exceed a defined distance cutoff.

    These criteria apply regardless of the order of the atoms in the pair, meaning a positively ionizable
    to negatively ionizable contact is considered equivalent to a negatively ionizable to positively ionizable contact.

    Attributes:
        ns (NeighborPairs): A NeighborPairs object containing the atom neighbor relationships in the system.
        distance (float): The maximum distance to consider for an ionic contact. See `lahuta.config.defaults.CONTACTS`
            for default values.

    ??? example "Example"
        ``` py
        luni = Luni(...)
        ns = luni.compute_neighbors()

        ionic = IonicContacts(ns)
        print(ionic.results)
        ```
    """

    distance = CONTACTS["ionic"]["distance"]

    def compute(self) -> NeighborPairs:
        """Compute ionic contacts based on the neighbor pairs.

        Returns:
            NeighborPairs: A NeighborPairs object containing only ionic contacts.

        """
        return F.ionic_neighbors(self.ns, self.distance)
