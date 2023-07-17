from abc import ABC, abstractmethod

from lahuta.core.neighbors import NeighborPairs


class ContactAnalysis(NeighborPairs):
    def __init__(self, ns):
        self.ns = ns
        super().__init__(ns.luni, ns.pairs, ns.distances)

        # self._results = self._compute(self.ns.luni.to("mda").atoms, self.pairs, self.distances)
        self._results = self._compute(
            self.ns.partner1, self.ns.partner2, self.distances
        )

    @abstractmethod
    def _compute(self, atoms, pairs, distances):
        ...
