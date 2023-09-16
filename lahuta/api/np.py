"""API for NeighborPairs objects."""
from typing import TypeVar, overload

from lahuta.core.labeled_neighbors import LabeledNeighborPairs
from lahuta.core.neighbors import NeighborPairs

__all__ = [
    "union",
    "intersection",
    "difference",
    "symmetric_difference",
]

T = TypeVar("T", NeighborPairs, LabeledNeighborPairs)


@overload
def two_element_union(a: NeighborPairs, b: NeighborPairs) -> NeighborPairs:
    ...


@overload
def two_element_union(a: LabeledNeighborPairs, b: LabeledNeighborPairs) -> LabeledNeighborPairs:
    ...


def two_element_union(a: T, b: T) -> T:
    """Compute the union of two neighbor type objects.

    Args:
        a (NeighborPairs | LabeledNeighborPairs): A NeighborPairs or LabeledNeighborPairs object.
        b (NeighborPairs | LabeledNeighborPairs): A NeighborPairs or LabeledNeighborPairs object.

    Returns:
        NeighborPairs | LabeledNeighborPairs: A NeighborPairs or LabeledNeighborPairs object.
    """
    return a + b


@overload
def two_element_intersection(a: NeighborPairs, b: NeighborPairs) -> NeighborPairs:
    ...


@overload
def two_element_intersection(a: LabeledNeighborPairs, b: LabeledNeighborPairs) -> LabeledNeighborPairs:
    ...


def two_element_intersection(a: T, b: T) -> T:
    """Compute the intersection of two neighbor type objects.

    Args:
        a (NeighborPairs | LabeledNeighborPairs): A NeighborPairs or LabeledNeighborPairs object.
        b (NeighborPairs | LabeledNeighborPairs): A NeighborPairs or LabeledNeighborPairs object.

    Returns:
        NeighborPairs | LabeledNeighborPairs: A NeighborPairs or LabeledNeighborPairs object.
    """
    return a & b


@overload
def two_element_difference(a: NeighborPairs, b: NeighborPairs) -> NeighborPairs:
    ...


@overload
def two_element_difference(a: LabeledNeighborPairs, b: LabeledNeighborPairs) -> LabeledNeighborPairs:
    ...


def two_element_difference(a: T, b: T) -> T:
    """Compute the difference of two neighbor type objects.

    Args:
        a (NeighborPairs | LabeledNeighborPairs): A NeighborPairs or LabeledNeighborPairs object.
        b (NeighborPairs | LabeledNeighborPairs): A NeighborPairs or LabeledNeighborPairs object.

    Returns:
        NeighborPairs | LabeledNeighborPairs: A NeighborPairs or LabeledNeighborPairs object.
    """
    return a - b


@overload
def two_element_symmetric_difference(a: NeighborPairs, b: NeighborPairs) -> NeighborPairs:
    ...


@overload
def two_element_symmetric_difference(a: LabeledNeighborPairs, b: LabeledNeighborPairs) -> LabeledNeighborPairs:
    ...


def two_element_symmetric_difference(a: T, b: T) -> T:
    """Compute the symmetric difference of two neighbor type objects.

    Args:
        a (NeighborPairs | LabeledNeighborPairs): A NeighborPairs or LabeledNeighborPairs object.
        b (NeighborPairs | LabeledNeighborPairs): A NeighborPairs or LabeledNeighborPairs object.

    Returns:
        NeighborPairs | LabeledNeighborPairs: A NeighborPairs or LabeledNeighborPairs object.
    """
    return a | b


@overload
def union(*args: NeighborPairs) -> NeighborPairs:
    ...


@overload
def union(*args: LabeledNeighborPairs) -> LabeledNeighborPairs:
    ...


def union(*args: T) -> T:
    """Compute the union of NeighborPairs or LabeledNeighborPairs objects.

    Args:
        *args (NeighborPairs | LabeledNeighborPairs): NeighborPairs or LabeledNeighborPairs objects.

    Returns:
        NeighborPairs | LabeledNeighborPairs: A NeighborPairs or LabeledNeighborPairs object.

    Examples:
        >>> union(ns1, ns2, ns3, ...)
        >>> ns_list = [ns1, ns2, ns3, ...]
        >>> union(*ns_list)
    """
    if len(args) == 0 or (len(args) == 1 and not args[0]):
        raise ValueError("No arguments provided.")

    result = args[0]
    for element in args[1:]:
        result += element

    return result


@overload
def intersection(*args: NeighborPairs) -> NeighborPairs:
    ...


@overload
def intersection(*args: LabeledNeighborPairs) -> LabeledNeighborPairs:
    ...


def intersection(*args: T) -> T:
    """Compute the intersection of NeighborPairs or LabeledNeighborPairs objects.

    Args:
        *args (NeighborPairs | LabeledNeighborPairs): NeighborPairs or LabeledNeighborPairs objects.

    Returns:
        NeighborPairs | LabeledNeighborPairs: A NeighborPairs or LabeledNeighborPairs object.

    Examples:
        >>> intersection(ns1, ns2, ns3, ...)
        >>> ns_list = [ns1, ns2, ns3, ...]
        >>> intersection(*ns_list)
    """
    if len(args) == 0 or (len(args) == 1 and not args[0]):
        raise ValueError("No arguments provided.")

    result = args[0]
    for element in args[1:]:
        result &= element

    return result


@overload
def difference(*args: NeighborPairs) -> NeighborPairs:
    ...


@overload
def difference(*args: LabeledNeighborPairs) -> LabeledNeighborPairs:
    ...


def difference(*args: T) -> T:
    """Compute the difference of NeighborPairs or LabeledNeighborPairs objects.

    Args:
        *args (NeighborPairs | LabeledNeighborPairs): NeighborPairs or LabeledNeighborPairs objects.

    Returns:
        NeighborPairs | LabeledNeighborPairs: A NeighborPairs or LabeledNeighborPairs object.

    Examples:
        >>> difference(ns1, ns2, ns3, ...)
        >>> ns_list = [ns1, ns2, ns3, ...]
        >>> difference(*ns_list)
    """
    if len(args) == 0 or (len(args) == 1 and not args[0]):
        raise ValueError("No arguments provided.")

    result = args[0]
    for element in args[1:]:
        result -= element

    return result


@overload
def symmetric_difference(*args: NeighborPairs) -> NeighborPairs:
    ...


@overload
def symmetric_difference(*args: LabeledNeighborPairs) -> LabeledNeighborPairs:
    ...


def symmetric_difference(*args: T) -> T:
    """Compute the symmetric difference of NeighborPairs or LabeledNeighborPairs objects.

    Args:
        *args (NeighborPairs | LabeledNeighborPairs): NeighborPairs or LabeledNeighborPairs objects.

    Returns:
        NeighborPairs | LabeledNeighborPairs: A NeighborPairs or LabeledNeighborPairs object.

    Examples:
        >>> symmetric_difference(ns1, ns2, ns3, ...)
        >>> ns_list = [ns1, ns2, ns3, ...]
        >>> symmetric_difference(*ns_list)
    """
    if len(args) == 0 or (len(args) == 1 and not args[0]):
        raise ValueError("No arguments provided.")

    result = args[0]
    for element in args[1:]:
        result |= element

    return result
