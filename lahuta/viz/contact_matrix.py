from typing import Literal, Tuple

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from .plotters import FullPlotter, MatchingIndicesPlotter


class ContactMap:
    def __init__(self, pairs: NDArray[np.int32], figsize: Tuple[int, int] = (10, 10)) -> None:
        self.pairs = pairs
        plt.figure(figsize=figsize)

    def plot(self, which: Literal['matching', 'full'] = 'matching', outline: bool = False, half_only: bool = False) -> None:
        if which == 'full':
            self.plot_full(outline, half_only)
        elif which == 'matching':
            self.plot_matching_indices()

    def plot_full(self, outline: bool = False, half_only: bool = False) -> None:
        plotter = FullPlotter(self.pairs)
        plotter.plot(outline, half_only)

    def plot_matching_indices(self) -> None:
        plotter = MatchingIndicesPlotter(self.pairs)
        plotter.plot()
