"""Contains the ContactMap class for plotting contact maps."""
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from .plotters import FullPlotter, MatchingIndicesPlotter


class ContactMap:
    """Plots a contact map visualizing the contacts between atoms.

    Args:
        pairs: A 2D array of shape (N, 2) where N is the number of pairs.
        figsize: The size of the figure to plot.
    """

    def __init__(self, pairs: NDArray[np.int32], figsize: tuple[int, int] = (10, 10)) -> None:
        self.pairs = pairs
        plt.figure(figsize=figsize)

    def plot(self, which: Literal["matching", "full"] = "matching", half_only: bool = False) -> None:
        """Plot the contact map.

        Args:
            which: Which contact map to plot. Either 'matching' or 'full'.
            half_only: Whether to plot only the upper half of the contact map.
        """
        if which == "full":
            self.plot_full(False, half_only)
        elif which == "matching":
            self.plot_matching_indices()

    def plot_full(self, outline: bool = False, half_only: bool = False) -> None:
        """Plot the full contact map.

        Args:
            outline: Whether to outline the contact map.
            half_only: Whether to plot only the upper half of the contact map.
        """
        plotter = FullPlotter(self.pairs)
        plotter.plot(outline, half_only)

    def plot_matching_indices(self) -> None:
        """Plot the contact map for only indices that are in contact."""
        plotter = MatchingIndicesPlotter(self.pairs)
        plotter.plot()
