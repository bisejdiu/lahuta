"""
Module: metal.py

This module defines a class for computing metal contacts using a class-based approach. 
The MetalicContacts class inherits from the base ContactAnalysis class and 
implements the `compute` method for metal contact computation.
                                       
Example:
    luni = Luni(...)
    ns = luni.compute_neighbors()

    metal_contacts = MetalicContacts(ns)
    print(metal_contacts.results)
"""

import lahuta.contacts as F
from lahuta.config.defaults import CONTACTS
from lahuta.core.neighbors import NeighborPairs

from .base import ContactAnalysis


class MetalicContacts(ContactAnalysis):
    """
    Handles the computation of metallic contacts in a molecular system.

    Metallic contacts are interactions involving metal atoms in a molecular system.
    This class, a derivative of the `ContactAnalysis` base class, overrides the `compute`
    method to provide functionality specifically for metallic contact computation.

    Metallic contacts are interactions between metal ions and atoms that can act as hydrogen bond acceptors.
    These contacts play significant roles in the structure and function of many proteins, especially metalloproteins.
    Metal ions can form coordination bonds with electron-rich atoms (like oxygen, nitrogen, or sulfur),
    contributing to the structural stability and sometimes the catalytic activity of these proteins.


    !!! tip "Definition"
        1. The contact involves a metal ion and an atom that is a hydrogen bond acceptor.
        2. The distance between the metal ion and the hydrogen bond acceptor does not exceed a predefined distance cutoff.

    Attributes:
        ns (NeighborPairs): A NeighborPairs object containing the atom neighbor relationships in the system.
        distance (float): The maximum distance to consider for a metallic contact. See `lahuta.config.defaults.CONTACTS`
            for default values.

    ??? example "Example"
        ``` py
        luni = Luni(...)
        ns = luni.compute_neighbors()

        metal_contacts = MetalicContacts(ns)
        print(metal_contacts.results)
        ```
    """

    distance = CONTACTS["metal"]["distance"]

    def compute(self) -> NeighborPairs:
        """Computes metalic contacts based on the neighbor pairs.

        Returns:
            NeighborPairs: A NeighborPairs object containing only metalic contacts.
        """
        return F.metalic_neighbors(self.ns, self.distance)
