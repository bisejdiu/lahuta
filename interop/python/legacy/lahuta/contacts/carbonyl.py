"""Defines a class for computing carbonyl contacts using a class-based approach.
The CarbonylContacts class inherits from the base ContactAnalysis class and 
implements the `compute` method for carbonyl contact computation.

Class:
    CarbonylContacts(ContactAnalysis): Computes carbonyl contacts.
                                       

Example:
    ``` py
    luni = Luni(...)
    ns = luni.compute_neighbors()

    cbyl = CarbonylContacts(ns)
    print(cbyl.results)
    ```

"""


import lahuta.contacts as F
from lahuta.config.defaults import CONTACTS
from lahuta.core.neighbors import NeighborPairs

from .base import ContactAnalysis


class CarbonylContacts(ContactAnalysis):
    """Handle the computation of carbonyl contacts in a molecular system.

    Carbonyl contacts involve the interaction between a carbonyl oxygen atom (O) and a carbonyl carbon atom (C)
    from a carbonyl functional group (C=O) in the context of protein-ligand structures or protein-protein structures.

    In a carbonyl group, the carbon atom has a double bond with the oxygen atom. This arrangement results in a polar
    bond with the oxygen atom carrying a partial negative charge and the carbon atom a partial positive charge.
    This polarity can lead to interactions with other polar or charged atoms.


    !!! tip "Definition"
        1. One atom is a carbonyl oxygen atom (O=C).
        2. The second atom is a carbonyl carbon atom (O=C).
        3. The distance between these two atoms does not exceed a defined distance cutoff.

    The directionality of the contact is not considered, meaning that an Oxygen to Carbon contact is equivalent
    to a Carbon to Oxygen contact.

    Attributes:
        ns (NeighborPairs): A NeighborPairs object containing the atom neighbor relationships in the system.
        distance (float): The maximum distance to consider for a carbonyl contact. See `lahuta.config.defaults.CONTACTS`
            for default values.

    ??? example "Example"
        ``` py
        luni = Luni(...)
        ns = luni.compute_neighbors()

        cbyl = CarbonylContacts(ns)
        print(cbyl.results)
        ```
    """

    distance = CONTACTS["carbonyl"]["distance"]

    def compute(self) -> NeighborPairs:
        """Compute carbonyl contacts based on the neighbor pairs.

        Returns:
            NeighborPairs: A NeighborPairs object containing only carbonyl contacts.
        """
        return F.carbonyl_neighbors(self.ns, self.distance)
