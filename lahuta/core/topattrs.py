"""
Module: topattrs.py

This module defines a class for dynamically generating new AtomAttr subclasses.

"""
from typing import Any, Optional, Type

import numpy as np
from MDAnalysis.core.topologyattrs import AtomAttr
from numpy.typing import NDArray


class AtomAttrClassHandler:
    """
    A class for dynamically generating new AtomAttr subclasses.

    This class facilitates the dynamic generation of new subclasses of AtomAttr, based
    on provided attribute and singular names. The attribute name informs the creation of
    the class attribute, while the singular name is used to derive the new class name.

    For instance, given an attribute name "atom_types" and a singular name "atom_type",
    a new subclass named "AtomTypes" with a class attribute "atom_type" will be created.

    Attributes:
        atomattr_class (AtomAttr subclass): The dynamically generated subclass of AtomAttr.
    """

    def __init__(self) -> None:
        self.atomattr_class: Optional[Type[AtomAttr]] = None

    @staticmethod
    def _gen_initial_values(n_atoms: int, *_: Any) -> NDArray[Any]:
        """
        Generate the initial values for the attribute.

        Args:
            n_atoms (int): The number of atoms.
            *_ (Any): Any other arguments are ignored.

        Returns:
            NDArray[Any]: A numpy array of zeros with a size equal to the number of atoms.
        """
        return np.zeros(n_atoms)

    def init_topattr(self, attrname: str, singular_name: str) -> None:
        """
        Generates a new AtomAttr subclass and adds it to the global namespace.

        Args:
            attrname (str): The name of the attribute.
            singular_name (str): The singular name of the attribute, used to generate the class name.
        """

        attr_dict = {
            "attrname": attrname,
            "singular": singular_name,
            "_gen_initial_values": self._gen_initial_values,
        }
        self.atomattr_class = type(attrname, (AtomAttr,), attr_dict)
        globals()[attrname] = self.atomattr_class


class VDWRadiiAtomAttr(AtomAttr):
    """
    A statically generated AtomAttr subclass for van der Waals radii.

    The reason for the static definition even though `AtomAttrClassHandler` can dynamically
    generate subclasses is that the van der Waals radii are used during contact calculations and
    parallelization is not possible with dynamically generated subclasses.

    Attributes:
        attrname (str): The name of the attribute.
        singular (str): The singular name of the attribute.
    """

    attrname = "vdw_radii"
    singular = "vdw_radii"

    @staticmethod
    def _gen_initial_values(n_atoms: int, *_: Any) -> NDArray[Any]:
        """
        Generate the initial values for the attribute.

        Args:
            n_atoms (int): The number of atoms.
            *_ (Any): Any other arguments are ignored.

        Returns:
            NDArray[Any]: A numpy array of zeros with a size equal to the number of atoms.
        """
        return np.zeros(n_atoms)
