"""Lahuta is a Python library for analyzing contacts in structural data and molecular dynamics trajectories."""

from lahuta.core import Luni, NeighborPairs
# from lahuta.lib import cLuni

from ._version import VERSION

__version__ = VERSION

__all__ = ["Luni", "NeighborPairs"]  # , "cLuni"]
__all__ += ["VERSION", "__version__"]

# Import mapping submodule
try:
    from . import mapping
    __all__.append("mapping")
except ImportError:
    pass
