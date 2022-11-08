"""
Placeholder for Obabel object.
"""

from openbabel import openbabel as ob


class OBMol:
    """A wrapper around openbabel.OBMol that adds some extra functionality."""

    def __init__(self, file_name):
        ob_conv = ob.OBConversion()
        ob_conv.SetInFormat("pdb")
        self.mol = ob.OBMol()
        ob_conv.ReadFile(self.mol, file_name)

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
