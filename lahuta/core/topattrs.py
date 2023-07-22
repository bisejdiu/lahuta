from typing import Any, Optional, Type

import numpy as np
from MDAnalysis.core.topologyattrs import AtomAttr
from numpy.typing import NDArray


class AtomAttrClassHandler:
    """A class for generating new AtomAttr subclasses.

    This class is used to generate new AtomAttr subclasses and add them to the
    global namespace. The new subclasses are generated based on the name of the
    attribute and the singular name of the attribute. The singular name is used
    to generate the name of the class, while the name of the attribute is used
    to generate the name of the class attribute.

    For example, if the attribute name is "atom_types" and the singular name is
    "atom_type", the new class will be named "AtomTypes" and the class attribute
    will be named "atom_type".

    Attributes:
        atomattr_class (AtomAttr subclass): The generated AtomAttr subclass.
    """

    def __init__(self) -> None:
        """Initialize the handler."""
        self.atomattr_class: Optional[Type[AtomAttr]] = None

    @staticmethod
    def _gen_initial_values(n_atoms: int, *_: Any) -> NDArray[Any]:
        """Generate the initial values for the attribute."""
        return np.zeros(n_atoms)

    def init_topattr(self, attrname: str, singular_name: str) -> None:
        """Generate the new AtomAttr subclass and add it to the global namespace.

        Args:
            attr_name (str): The name of the attribute.
            singular_name (str): The singular name of the attribute.

        Returns:
            AtomAttr subclass: The created class is added to the global namespace.
        """

        attr_dict = {
            "attrname": attrname,
            "singular": singular_name,
            "_gen_initial_values": self._gen_initial_values,
        }
        self.atomattr_class = type(attrname, (AtomAttr,), attr_dict)
        globals()[attrname] = self.atomattr_class
