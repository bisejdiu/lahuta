"""
Module: mda_commands.py

This module provides typed wrappers for certain MDAnalysis commands, which don't natively support typing. 
The primary objective of this module is to enhance static type checking, providing an additional 
layer of error detection and contributing to cleaner code.

This module currently focuses on commands related to computing distances, specifically the `capped_distance` 
command. The Protocol-based DistanceType class is used to ensure typing for the command, and the 
CappedDistance class is a wrapper implementing this command with appropriate type hints.

As MDAnalysis evolves to incorporate typing natively, modules like this will be phased out, simplifying the codebase.

Classes:
    DistanceType(Protocol): Provides a typing interface for distance-related commands in MDAnalysis.
    CappedDistance: Wraps the `capped_distance` command from MDAnalysis with appropriate type hints.

Example:
    capped_distance_command = ...
    cd = CappedDistance(capped_distance_command)
    reference = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    configuration = np.array([[1.0, 1.0, 1.0]], dtype=np.float32)
    pairs, distances = cd.capped_distance(reference, configuration, 1.5)
    print(pairs, distances)
"""

from typing import Any, Protocol, Tuple

import numpy as np
from numpy.typing import NDArray


# pylint: disable=missing-function-docstring
# pylint: disable=too-few-public-methods
class DistanceType(Protocol):
    """A typing interface for distance-related commands in MDAnalysis."""

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


# pylint: disable=too-many-arguments
# pylint: disable=too-few-public-methods
class CappedDistance:
    """A wrapper for the `capped_distance` command from MDAnalysis with appropriate type hints.

    Args:
        capped_distance: The `capped_distance` command from MDAnalysis.

    """

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
