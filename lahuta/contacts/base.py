"""Defines base classes and protocols for implementing atomic contact computation using a class-based approach.
It provides two main classes: ComputeProtocol and ComputeElementwiseProtocol, and a base class ContactAnalysis which 
defines the core logic for running contact computation methods.

Protocols:
    ComputeProtocol: A runtime-checkable protocol which requires classes to implement a `compute` method. 

    ComputeElementwiseProtocol: A runtime-checkable protocol which requires classes to implement 
    a `compute_elementwise` method.

Class:
    ContactAnalysis: A base class that provides the structure for running contact computation methods. 
                     Contact computation classes should inherit from this class and implement methods 
                     as required by the ComputeProtocol or ComputeElementwiseProtocol.

To compute contacts using the class-based approach, create a new class inheriting from ContactAnalysis and 
implement the `compute` method or `compute_elementwise` method as per requirements.

Example:
    class CustomContacts(ContactAnalysis):
        distance: float = 4.5
        def compute(self) -> NeighborPairs:
            "Compute contacts based on a custom distance."
            return self.ns.distance_filter(self.distance)

    class ElementwiseCustomContacts(ContactAnalysis):
        distance: float = 4.5
        def compute_elementwise(self, atom1, atom2, distance) -> bool:
            "Compute contacts by checking each atom pair one-by-one."
            return distance <= self.distance

Notes:
    Implementing the `compute` method enables very fast computation of contacts as it works on the entire 
    NeighborPairs instance. However, it requires dealing with NeighborPairs arithmetic operations and assumes 
    familiarity with the overall system and neighbor pairs.

    On the other hand, implementing `compute_elementwise` method can be slower as it checks each atom pair 
    individually, but it is easier to work with, especially when custom conditions involving specific atom 
    attributes need to be checked for contact computation.

"""

from typing import Any, Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray

from lahuta.api import union
from lahuta.core.neighbors import NeighborPairs
from lahuta.lahuta_types.mdanalysis import AtomGroupType

T = Any


@runtime_checkable
class ComputeProtocol(Protocol):
    """A runtime-checkable protocol which requires classes to implement a `compute` method."""

    def compute(self) -> NeighborPairs:
        """Compute contacts based on the neighbor pairs."""

    def _conclude(self) -> T:
        """Conclude the computation and store the results."""

@runtime_checkable
class ComputeElementwiseProtocol(Protocol):
    """A runtime-checkable protocol which requires classes to implement a `compute_elementwise` method."""

    def compute_elementwise(
        self,
        atoms: AtomGroupType,
        pair: NDArray[np.int32],
        distance: NDArray[np.float_],
    ) -> T | None:
        """Compute contacts based on the neighbor pairs."""

    def _conclude(self) -> T:
        """Conclude the computation and store the results."""

class ContactAnalysis:
    """A base class that provides the structure for running contact computation methods.

    Contact computation classes should inherit from this class and implement methods as required by the
    ComputeProtocol or ComputeElementwiseProtocol.

    Attributes:
        ns (NeighborPairs): A NeighborPairs object containing the atom neighbor relationships in the system.
        results (Any): The results of the contact computation.
    """

    def __init__(self, ns: NeighborPairs):
        self.ns = ns
        self._results: Any = None

        self.run()
        self._conclude()

    def run(self) -> None:
        """Run the contact computation methods."""
        self.run_methods()

    def run_methods(self) -> None:
        """Run the compute or compute_elementwise methods."""
        if isinstance(self, ComputeProtocol):
            self._results = self.compute()
        elif isinstance(self, ComputeElementwiseProtocol):
            self._results = []
            p1_atoms, p2_atoms = (
                self.ns.partner1,
                self.ns.partner2,
            )
            for atom1, atom2, distance in zip(p1_atoms, p2_atoms, self.ns.distances, strict=True):
                result = self.compute_elementwise(atom1, atom2, distance)
                if result:
                    self._results.append(result)
        else:
            raise NotImplementedError("Object must implement either compute or compute_elementwise methods.")

    def _conclude(self) -> None:
        """Conclude the computation and store the results."""
        if isinstance(self, ComputeProtocol):
            self.results = self._results
        elif isinstance(self, ComputeElementwiseProtocol):
            self.results = union(*self._results) if self._results else NeighborPairs(self.ns.luni)
        else:
            raise NotImplementedError("Object must implement either compute or compute_elementwise methods.")