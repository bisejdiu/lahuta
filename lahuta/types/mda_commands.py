from typing import Any, Protocol, Tuple

import numpy as np
from numpy.typing import NDArray


class DistanceType(Protocol):
    def capped_distance(
        self,
        reference: NDArray[np.float32],
        configuration: NDArray[np.float32],
        max_cutoff: float,
        min_cutoff: None = None,
        box: None = None,
        method: None = None,
        return_distances: bool = True,
    ) -> Any:
        ...


class CappedDistance:
    def __init__(self, capped_distance: DistanceType):
        self.func = capped_distance

    def capped_distance(
        self,
        reference: NDArray[np.float32],
        configuration: NDArray[np.float32],
        max_cutoff: float,
        min_cutoff: None = None,
        box: None = None,
        method: None = None,
        return_distances: bool = True,
    ) -> Tuple[NDArray[np.int32], NDArray[np.float32]]:
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
