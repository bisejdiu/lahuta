from abc import ABC, abstractmethod

from openbabel import openbabel as ob


class FileLoader(ABC):
    def __init__(self, file_name: str, in_format: str = None):
        self.file_name = file_name
        self.in_format = in_format
        self.mol = None

    @abstractmethod
    def load(self, *args):
        raise NotImplementedError("Subclasses must implement this method")

    def _load_obabel(self):
        if self.in_format is None:
            try:
                self.in_format = self.file_name.split(".")[-1]
            except AttributeError as exc:
                raise AttributeError("You must either specify the input format using the in_format "
                                     "argument or pass a file name with a file extension.") from exc

        ob_conv = ob.OBConversion()
        ob_conv.SetInFormat(self.in_format)
        mol = ob.OBMol()
        ob_conv.ReadFile(mol, self.file_name)
        self.mol = mol

