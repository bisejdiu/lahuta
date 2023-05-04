"""
Placeholder for Obabel object.
"""

from openbabel import openbabel as ob


class OBMol:
    """A wrapper around openbabel.OBMol that adds some extra functionality."""

    def __init__(self, file_name: str, in_format: str = None):
        """A wrapper around openbabel.OBMol that adds some extra functionality."""
        if in_format is None:
            try:
                in_format = file_name.split(".")[-1]
            except AttributeError as exc:
                raise AttributeError(
                    "You must specify the input format if you are passing a string."
                ) from exc

        log_level = ob.cvar.obErrorLog.GetOutputLevel()
        ob.cvar.obErrorLog.SetOutputLevel(0)

        ob_conv = ob.OBConversion()
        ob_conv.SetInFormat(in_format)
        self.mol = ob.OBMol()
        # print ('Reading file:', file_name)
        ob_conv.ReadFile(self.mol, file_name)
        ob.cvar.obErrorLog.SetOutputLevel(log_level)

    def __getattr__(self, attr):
        return getattr(self.mol, attr)

    def __len__(self):
        return self.mol.NumAtoms()

    def __iter__(self):
        return ob.OBMolAtomIter(self.mol)

    def __str__(self):
        return f"OBMol with {len(self)} atoms"

    def __repr__(self):
        return self.__str__()
