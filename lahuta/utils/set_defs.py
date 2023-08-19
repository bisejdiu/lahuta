"""Set definitions for theorems.

This module defines theorems that can be checked for labeled neighbor pairs.
"""
import numpy as np
from numpy.typing import NDArray

from lahuta.core.labeled_neighbors import LabeledNeighborPairs


def check_union(s1: LabeledNeighborPairs, s2: LabeledNeighborPairs) -> bool:
    """Check Union theorem.

    Definition: |s1 U s2| = |s1| + |s2| - |s1 ∩ s2|

    Args:
        s1 (LabeledNeighborPairs): The first object to check.
        s2 (LabeledNeighborPairs): The second object to check.

    Returns:
        bool: True if the theorem holds, False otherwise.
    """
    return len(s1 + s2) == len(s1) + len(s2) - len(s1 & s2)


def check_intersection(s1: LabeledNeighborPairs, s2: LabeledNeighborPairs) -> bool:
    """Check Intersection theorem.

    Definition: |s1 ∩ s2| = |s1| + |s2| - |s1 U s2|

    Args:
        s1 (LabeledNeighborPairs): The first object to check.
        s2 (LabeledNeighborPairs): The second object to check.

    Returns:
        bool: True if the theorem holds, False otherwise.
    """
    return len(s1 & s2) == len(s1) + len(s2) - len(s1 + s2)


def check_difference(s1: LabeledNeighborPairs, s2: LabeledNeighborPairs) -> bool:
    """Check Difference theorem.

    Definition: |s1 - s2| = |s1| - |s1 ∩ s2|

    Args:
        s1 (LabeledNeighborPairs): The first object to check.
        s2 (LabeledNeighborPairs): The second object to check.

    Returns:
        bool: True if the theorem holds, False otherwise.
    """
    return len(s1 - s2) == len(s1) - len(s1 & s2)


def check_symmetric_difference(s1: LabeledNeighborPairs, s2: LabeledNeighborPairs) -> bool:
    """Check Symmetric Difference theorem.

    Definition: |s1 △ s2| = |s1| + |s2| - 2 * |s1 ∩ s2|

    Args:
        s1 (LabeledNeighborPairs): The first object to check.
        s2 (LabeledNeighborPairs): The second object to check.

    Returns:
        bool: True if the theorem holds, False otherwise.
    """
    return len(s1 | s2) == len(s1) + len(s2) - 2 * len(s1 & s2)


def check_subset(s1: LabeledNeighborPairs) -> bool:
    """Check Subset theorem.

    Definition: |s1| ≤ |s2| ⇔ s1 ⊆ s2

    Args:
        s1 (LabeledNeighborPairs): The first object to check.
        s2 (LabeledNeighborPairs): The second object to check.

    Returns:
        bool: True if the theorem holds, False otherwise.
    """
    m1_s = s1.create_new(random_subset(s1))
    return (len(m1_s) <= len(s1)) == (m1_s <= s1)


def check_superset(s1: LabeledNeighborPairs) -> bool:
    """Check Superset theorem.

    Definition: |s1| ≥ |s2| ⇔ s1 ⊇ s2

    Args:
        s1 (LabeledNeighborPairs): The first object to check.
        s2 (LabeledNeighborPairs): The second object to check.

    Returns:
        bool: True if the theorem holds, False otherwise.
    """
    m1_s = s1.create_new(random_subset(s1))
    return (len(m1_s) >= len(s1)) == (m1_s >= s1)


def check_proper_subset(s1: LabeledNeighborPairs) -> bool:
    """Check Proper Subset theorem.

    Definition: |s1| < |s2| ⇔ s1 ⊂ s2

    Args:
        s1 (LabeledNeighborPairs): The first object to check.
        s2 (LabeledNeighborPairs): The second object to check.

    Returns:
        bool: True if the theorem holds, False otherwise.
    """
    m1_s = s1.create_new(random_subset(s1))
    return (len(m1_s) < len(s1)) == (m1_s < s1)


def check_proper_superset(s1: LabeledNeighborPairs) -> bool:
    """Check Proper Superset theorem.

    Definition: |s1| > |s2| ⇔ s1 ⊃ s2

    Args:
        s1 (LabeledNeighborPairs): The first object to check.
        s2 (LabeledNeighborPairs): The second object to check.

    Returns:
        bool: True if the theorem holds, False otherwise.
    """
    m1_s = s1.create_new(random_subset(s1))
    return (len(m1_s) > len(s1)) == (m1_s > s1)


def check_disjoint(s1: LabeledNeighborPairs, s2: LabeledNeighborPairs) -> bool:
    """Check Disjoint theorem.

    Definition: |s1 ∩ s2| = 0 ⇔ s1 ∩ s2 = ∅

    Args:
        s1 (LabeledNeighborPairs): The first object to check.
        s2 (LabeledNeighborPairs): The second object to check.

    Returns:
        bool: True if the theorem holds, False otherwise.
    """
    # pylint: disable=import-outside-toplevel
    from lahuta.core.builder import LabeledNeighborPairsBuilder

    empty = s1.create_new(np.empty((0, 2), dtype=LabeledNeighborPairsBuilder.DTYPE))
    m1_subset = s1.create_new(random_subset(s1))
    m1_altered = m1_subset - s2
    return empty.isdisjoint(s2) == (len(empty & s2) == 0) and len(m1_altered & s2) == 0


def random_subset(s: LabeledNeighborPairs, size_ratio: float = 0.1) -> NDArray[np.void]:
    """Return a random subset of the given object.

    Args:
        s: The object from which the subset is to be taken.
        size_ratio (float, optional): The ratio of the size of the subset to the original object. Defaults to 0.1.

    Returns:
        The randomly selected subset of the given object.
    """
    size = s.pairs.shape[0]
    generator = np.random.default_rng()
    subset_indices = generator.choice(np.arange(size), size=int(size * size_ratio), replace=False)
    return s.pairs[subset_indices]  # type: ignore
