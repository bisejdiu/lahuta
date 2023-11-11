"""Lahuta is a Python library for analyzing contacts in structural data and molecular dynamics trajectories."""

from lahuta.core import Luni, NeighborPairs

from ._version import VERSION

__version__ = VERSION

__all__ = ["Luni", "NeighborPairs"]
__all__ += ["VERSION", "__version__"]
