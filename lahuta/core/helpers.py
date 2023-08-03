"""
This module contains helper functions for the lahuta.core package.

Functions:
    get_class_methods: Retrieves all the methods of the specified class.
    get_class_properties: Retrieves all the properties of the specified class.
    get_class_attributes: Retrieves all the attributes (not methods or properties) 
                        of the specified class.

"""
from typing import TYPE_CHECKING, List

if TYPE_CHECKING:
    from lahuta.core.neighbors import NeighborPairs


def get_class_methods(cls: "NeighborPairs") -> List[str]:
    """
    Retrieves all the methods of the specified class.

    Args:
        cls (NeighborPairs): The class to inspect.

    Returns:
        List[str]: A list of strings containing the names of all the methods of the specified class.
    """
    return [attr for attr in dir(cls) if callable(getattr(cls, attr))]


def get_class_properties(cls: "NeighborPairs") -> List[str]:
    """
    Retrieves all the properties of the specified class.

    Args:
        cls (NeighborPairs): The class to inspect.

    Returns:
        List[str]: A list of strings containing the names of all the properties of the specified class.
    """
    return [attr for attr in dir(cls) if isinstance(getattr(cls, attr), property)]


def get_class_attributes(cls: "NeighborPairs") -> List[str]:
    """
    Retrieves all the attributes (not methods or properties) of the specified class.

    Args:
        cls (NeighborPairs): The class to inspect.

    Returns:
        List[str]: A list of strings containing the names of all the attributes (not methods or properties)
        of the specified class.
    """
    return [attr for attr in dir(cls) if not attr.startswith("__") and not callable(getattr(cls, attr))]
