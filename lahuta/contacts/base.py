from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class ComputeProtocol(Protocol):
    def compute(self) -> Any:
        ...


@runtime_checkable
class ComputeElementwiseProtocol(Protocol):
    def compute_elementwise(self, atoms, pair, distance) -> Any:
        ...


class ContactAnalysis:
    def __init__(self, ns):
        self.ns = ns
        self.partner1_atoms = ns.atoms
        self.results: Any = None

        self.run()

    def run(self):
        self.run_methods()

    def run_methods(self):
        if isinstance(self, ComputeProtocol):
            self.results = self.compute()
        elif isinstance(self, ComputeElementwiseProtocol):
            self.results = []
            p1_atoms, p2_atoms = self.ns.partner1, self.ns.partner2
            for atom1, atom2, distance in zip(p1_atoms, p2_atoms, self.ns.distances):
                self.results.append(self.compute_elementwise(atom1, atom2, distance))
        else:
            raise NotImplementedError(
                "Object must implement either compute or compute_elementwise methods."
            )
