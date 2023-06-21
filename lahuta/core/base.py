from abc import ABC, abstractmethod
from typing import Optional

from openbabel import openbabel as ob


class FileLoader(ABC):
    """
    Base class for file loaders.

    Parameters
    ----------
    file_name : str
        The name of the file to load.

    in_format : str, optional
        The format of the input file. If not specified, the format will be
        inferred from the file extension.
    """

    def __init__(self, file_name: str, in_format: Optional[str] = None):
        self.file_name = file_name
        self.in_format = (
            in_format if in_format is not None else self.file_name.split(".")[-1]
        )
        self.mol = None

    @abstractmethod
    def load(self, *args):
        """Load the file and return the universe"""
        raise NotImplementedError("Subclasses must implement this method")

    def _load_obabel(self):
        ob_conv = ob.OBConversion()
        ob_conv.SetInFormat(self.in_format)
        mol = ob.OBMol()
        ob_conv.ReadFile(mol, self.file_name)
        self.mol = mol
