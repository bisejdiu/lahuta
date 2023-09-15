"""Type hints for gemmi objects."""
from typing import Protocol


class Structure(Protocol):
    
    def assign_serial_numbers(self) -> None:
        ...

    @property
    def cell(self) -> None: # actually it's a UnitCell instance that's returned
        ... 

    def __getitem__(self, index: int | slice) -> "Model":
        ...


class Model(Protocol):
    pass 

class SearchResults(Protocol):
    dist: float
    image_idx: int
    partner1: "GemmiCRA"
    partner2: "GemmiCRA"


class GemmiCRA(Protocol):
    atom: "GemmiAtom"

class GemmiAtom(Protocol):
    serial: int

    

