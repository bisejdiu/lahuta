"""Defines a class for computing covalent contacts using a class-based approach.
The CovalentContacts class inherits from the base ContactAnalysis class and 
implements the `compute` method for covalent contact computation.

Note:
    This contact is slightly slower because it requires a loop over all atoms in the molecule 
    to check for pairs that are covalently bonded. 
                               

Example:
    ``` py
    universe = Universe(...)
    ns = universe.compute_neighbors()

    cov = CovalentContacts(ns)
    print(cov.results)
    ```

"""


import lahuta.contacts as F
from lahuta.core.neighbors import NeighborPairs

from .base import ContactAnalysis


class CovalentContacts(ContactAnalysis):
    """Handle the computation of covalent contacts in a molecular system.

    Covalent contacts are interactions based on covalent bonds between atoms in a molecular system.
    This class, a derivative of the `ContactAnalysis` base class, overrides the `compute` method
    to provide functionality specifically for covalent contact computation.

    Covalent contacts refer to the interactions between atoms that share an electron pair, forming a covalent bond.
    We use the OpenBabel library to identify covalent bonds in the structure.

    !!! tip "Definition"
        Two atoms are considered to form a covalent contact if they are covalently bonded according to the molecular
        structure information obtained from OpenBabel.

    Attributes:
        ns (NeighborPairs): A NeighborPairs object containing the atom neighbor relationships in the system.


    ??? example "Example"
        ``` py
        luni = Luni(...)
        ns = luni.compute_neighbors()

        cov = CovalentContacts(ns)
        print(cov.results)
        ```

    """

    def compute(self) -> NeighborPairs:
        """Compute covalent contacts based on the neighbor pairs.

        Returns:
            NeighborPairs: A NeighborPairs object containing only covalent contacts.
        """
        return F.covalent_neighbors(self.ns)
