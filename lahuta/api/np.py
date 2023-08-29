"""API for NeighborPairs objects."""
from lahuta.core.neighbors import NeighborPairs

__all__ = [
    "union",
    "intersection",
    "difference",
    "symmetric_difference",
]


def union(a: NeighborPairs, b: NeighborPairs) -> NeighborPairs:
    """Compute the union of two NeighborPairs objects."""
    return a + b


def intersection(a: NeighborPairs, b: NeighborPairs) -> NeighborPairs:
    """Compute the intersection of two NeighborPairs objects."""
    return a & b


def difference(a: NeighborPairs, b: NeighborPairs) -> NeighborPairs:
    """Compute the difference of two NeighborPairs objects."""
    return a - b


def symmetric_difference(a: NeighborPairs, b: NeighborPairs) -> NeighborPairs:
    """Compute the symmetric difference of two NeighborPairs objects."""
    return a | b
