"""
Radius-based neighbor search.
"""
from __future__ import annotations
import numpy
import typing
__all__: list[str] = ['radius_neighbors', 'radius_neighbors_flat']
def radius_neighbors(X: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]], radius: float, Y: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]] | None = None, return_distance: bool = False, sort_results: bool = False) -> typing.Any:
    """
    Neighbors within radius. Returns lists per sample.
    """
def radius_neighbors_flat(X: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]], radius: float, Y: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]] | None = None, return_distance: bool = False, sort_results: bool = False) -> tuple:
    """
    Neighbors within radius. Returns CSR-like flat arrays: (indices, indptr) or (distances, indices, indptr).
    """
