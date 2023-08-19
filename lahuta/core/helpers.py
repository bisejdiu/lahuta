"""Helper functions to retrieve class methods and/or attributes. It is intended to be
used by and for NeighborPairs, even though the functions are generic.                
"""
from typing import TYPE_CHECKING, List

if TYPE_CHECKING:
    from lahuta.core.neighbors import NeighborPairs


def get_class_methods(cls: "NeighborPairs") -> List[str]:
    """Retrieve all the methods of the specified class.

    Args:
        cls (NeighborPairs): The class to inspect.

    Returns:
        List[str]: A list of strings containing the names of all the methods of the specified class.
    """
    return [attr for attr in dir(cls) if callable(getattr(cls, attr))]


def get_class_properties(cls: "NeighborPairs") -> List[str]:
    """Retrieve all the properties of the specified class.

    Args:
        cls (NeighborPairs): The class to inspect.

    Returns:
        List[str]: A list of strings containing the names of all the properties of the specified class.
    """
    return [attr for attr in dir(cls) if isinstance(getattr(cls, attr), property)]


def get_class_attributes(cls: "NeighborPairs") -> List[str]:
    """Retrieve all the attributes (not methods or properties) of the specified class.

    Args:
        cls (NeighborPairs): The class to inspect.

    Returns:
        List[str]: A list of strings containing the names of all the attributes (not methods or properties)
        of the specified class.
    """
    return [attr for attr in dir(cls) if not attr.startswith("__") and not callable(getattr(cls, attr))]
