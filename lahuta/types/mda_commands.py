from typing import Any, Callable, List, Protocol, Tuple

import numpy as np
from numpy.typing import NDArray


class MDADistanceType(Protocol):
    def capped_distance(
        self,
        reference: NDArray[np.float_],
        configuration: NDArray[np.float_],
        max_cutoff: float,
        min_cutoff: None = None,
        box: None = None,
        method: None = None,
        return_distances: bool = True,
    ) -> Any:
        ...


class MDACappedDistance:
    def __init__(self, capped_distance: MDADistanceType):
        self.func = capped_distance

    def capped_distance(
        self,
        reference: NDArray[np.float_],
        configuration: NDArray[np.float_],
        max_cutoff: float,
        min_cutoff: None = None,
        box: None = None,
        method: None = None,
        return_distances: bool = True,
    ) -> Tuple[NDArray[np.int_], NDArray[np.float_]]:
        pairs, distances = self.func.capped_distance(
            reference,
            configuration,
            max_cutoff,
            min_cutoff,
            box,
            method,
            return_distances,
        )
        return pairs, distances
