"""
Module: hydrophobic.py

This module defines a class for computing hydrophobic contacts using a class-based approach. 
The HydrophobicContacts class inherits from the base ContactAnalysis class and 
implements the `compute` method for hydrophobic contact computation.

Class:
    HydrophobicContacts(ContactAnalysis): Computes hydrophobic contacts.
                                       
Example:
    luni = Luni(...)
    ns = luni.compute_neighbors()

    hbc = HydrophobicContacts(ns)
    print(hbc.results)
"""

import lahuta.contacts as F
from lahuta.config.defaults import CONTACTS
from lahuta.core.neighbors import NeighborPairs

from .base import ContactAnalysis


class HydrophobicContacts(ContactAnalysis):
    """
    Handles the computation of hydrophobic contacts in a molecular system.

    Hydrophobic contacts are interactions between hydrophobic atoms in a molecular system.
    This class, a derivative of the `ContactAnalysis` base class, overrides the `compute`
    method to provide functionality specifically for hydrophobic contact computation.

    Hydrophobic contacts are interactions between hydrophobic (non-polar) atoms.
    Hydrophobic interactions occur due to the tendency of hydrophobic molecules to aggregate together
    in an aqueous environment, minimizing their exposure to water molecules.

    !!! tip "Definition"
        1. Both atoms must be hydrophobic.
        2. The distance between these two atoms does not exceed a defined distance cutoff.

    Attributes:
        ns (NeighborPairs): A NeighborPairs object containing the atom neighbor relationships in the system.
        distance (float): The maximum distance to consider for a hydrophobic contact. See `lahuta.config.defaults.CONTACTS`
            for default values.

    ??? example "Example"
        ``` py
        luni = Luni(...)
        ns = luni.compute_neighbors()

        hbc = HydrophobicContacts(ns)
        print(hbc.results)
        ```
    """

    distance = CONTACTS["hydrophobic"]["distance"]

    def compute(self) -> NeighborPairs:
        """Computes hydrophobic contacts based on the neighbor pairs.

        Returns:
            NeighborPairs: A NeighborPairs object containing only hydrophobic contacts.
        """
        return F.hydrophobic_neighbors(self.ns, self.distance)
