"""Defines a class for computing vanderwaals contacts using a class-based approach.
The VanDerWaalsContacts class inherits from the base ContactAnalysis class and 
implements the `compute` method for vanderwaals contact computation.

Class:
    VanDerWaalsContacts(ContactAnalysis): Computes vanderwaals contacts.
                                       

Example:
    universe = Universe(...)
    ns = universe.compute_neighbors()

    vdw = VanDerWaalsContacts(ns)
    print(vdw.results)

"""

import lahuta.contacts as F
from lahuta.core.neighbors import NeighborPairs

from .base import ContactAnalysis


class VanDerWaalsContacts(ContactAnalysis):
    """Handle the computation of Van der Waals (VdW) contacts in a molecular system.

    Van der Waals (VdW) contacts are determined based on the interactions between atoms that come
    within their combined van der Waals radii, increased by a compensation factor.
    Van der Waals interactions occur due to induced polarization of atoms and can play a critical role
    in the stability of biological structures and in molecular recognition processes.

    !!! tip "Definition"
        1. The distance between two atoms does not exceed the sum of their van der Waals radii,
        increased by a defined compensation factor.
        2. Optionally, if the distance between two atoms is less than the sum of their van der Waals radii
        (defined as a clash), the contact can be excluded.

    Attributes:
        ns (NeighborPairs): A NeighborPairs object containing the atom neighbor relationships in the system.
        vdw_comp_factor (float): The factor by which Van der Waals radii are multiplied during contact computation.
            This factor is used in the F.vdw_neighbors function to scale the VdW radii. Default value is 0.1.
        remove_clashes (bool): Flag indicating whether clashes (interactions with distances smaller than
            the combined Van der Waals radii of two atoms) should be removed. Default value is True.

    """

    vdw_comp_factor = 0.1
    remove_clashes = True

    def compute(self) -> NeighborPairs:
        """Compute Van der Waals contacts based on the neighbor pairs.

        Returns:
            NeighborPairs: A NeighborPairs object containing only Van der Waals contacts.
        """
        return F.vdw_neighbors(self.ns, self.vdw_comp_factor, self.remove_clashes)
