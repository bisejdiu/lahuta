"""Classes for plotting contact maps."""
import matplotlib.pyplot as plt
import numpy as np

from .base import BasePlotter


class FullPlotter(BasePlotter):
    """Plot the full contact map.

    Args:
        pairs: A list of pairs of indices that should be plotted.
    """

    def plot(self, outline: bool = False, half_only: bool = False) -> None:
        """Plot the full contact map.

        Args:
            outline: Whether to outline the contact map.
            half_only: Whether to plot only the upper half of the contact map.
        """
        self.half_only = half_only

        min_idx = np.min(self.pairs)
        max_idx = np.max(self.pairs)
        shape = (max_idx - min_idx + 1, max_idx - min_idx + 1)
        offset_pairs = [(x - min_idx, y - min_idx) for x, y in self.pairs]
        contact_map = self._initialize_map(shape, offset_pairs)

        plt.imshow(contact_map, cmap=self.binary_cmap, interpolation="none", origin="lower")
        if outline:
            self._add_outlines(min_idx, min_idx)

        plt.xlabel("Atom Indices")
        plt.ylabel("Atom Indices")
        plt.show()

    def _add_outlines(self, x_offset: np.int32, y_offset: np.int32) -> None:
        unique_x = np.unique(self.pairs[:, 0]) - x_offset
        unique_y = np.unique(self.pairs[:, 1]) - y_offset
        for x in unique_x:
            plt.axhline(y=x, color="grey", linestyle="--", lw=0.3)
        for y in unique_y:
            plt.axvline(x=y, color="grey", linestyle="--", lw=0.3)


class MatchingIndicesPlotter(BasePlotter):
    """Plot the contact map for only indices that are in contact.

    Args:
        pairs: A list of pairs of indices that should be plotted.
    """

    def plot(self) -> None:
        """Plot the contact map for only indices that are in contact."""
        self.half_only = True
        x_indices = np.unique(self.pairs[:, 0])
        y_indices = np.unique(self.pairs[:, 1])
        shape = (len(x_indices), len(y_indices))

        x_mapped = np.searchsorted(x_indices, self.pairs[:, 0])
        y_mapped = np.searchsorted(y_indices, self.pairs[:, 1])

        mapped_pairs = list(zip(x_mapped, y_mapped))

        contact_map = self._initialize_map(shape, mapped_pairs)
        plt.imshow(contact_map, cmap=self.binary_cmap, interpolation="none", origin="lower")

        n_points = min(x_indices.shape[0], y_indices.shape[0])
        if n_points < 10:
            plt.xticks(ticks=np.arange(len(y_indices)), labels=y_indices)
            plt.yticks(ticks=np.arange(len(x_indices)), labels=x_indices)
        else:
            x_tick_positions, x_tick_labels = self.get_ticks_and_labels(x_indices, n_points)
            y_tick_positions, y_tick_labels = self.get_ticks_and_labels(y_indices, n_points)

            plt.xticks(ticks=x_tick_positions, labels=x_tick_labels)
            plt.yticks(ticks=y_tick_positions, labels=y_tick_labels)

        plt.xlabel("Atom Indices")
        plt.ylabel("Atom Indices")
        plt.show()
