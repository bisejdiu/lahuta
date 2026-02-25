"""
Distance metrics for 3D coordinates.
"""
from __future__ import annotations
import typing
__all__: list[str] = ['cdist', 'pairwise_distances', 'pdist']
def cdist(XA: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]], XB: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]], squared: bool = False) -> numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]:
    """
    Compute distances between two sets of 3D points.
    """
def pairwise_distances(X: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]], Y: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]] | None = None, squared: bool = False) -> numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]:
    """
    Compute pairwise Euclidean distances.
    """
def pdist(X: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]], squared: bool = False) -> numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]:
    """
    Pairwise distances for one set in condensed form.
    """
