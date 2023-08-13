from typing import List, Tuple, Union

import numpy as np
from matplotlib.colors import ListedColormap
from numpy.typing import NDArray

AnyInt = Union[int, np.int32]

class BasePlotter:
    def __init__(self, pairs: NDArray[np.int32]) -> None:
        self.pairs = pairs
        self.half_only = False
        self.binary_cmap = ListedColormap(['white', 'black']) # type: ignore

    def plot(self) -> None:
        raise NotImplementedError()

    def _initialize_map(self, shape: Tuple[AnyInt, AnyInt], pairs: List[Tuple[np.int32, np.int32]]) -> NDArray[np.int32]:
        contact_map = np.zeros(shape, dtype=int)
        for x, y in pairs:
            contact_map[x, y] = 1
            if not self.half_only:
                contact_map[y, x] = 1
        return contact_map

    @staticmethod
    def get_ticks_and_labels(indices: NDArray[np.int32], n_points: int) -> Tuple[NDArray[np.int32], NDArray[np.int32]]:
        """
        Returns the positions and labels of the ticks for the given indices.
        The number of ticks is limited to 10.
        """
        n_points = min(10, n_points)
        if n_points == 0: 
            return np.array([]), np.array([])
        tick_positions = np.linspace(0, len(indices)-1, n_points).astype(int)
        tick_labels = indices[tick_positions]
        return tick_positions, tick_labels

