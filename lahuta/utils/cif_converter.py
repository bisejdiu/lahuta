"""
Converts CIF files to PDB files. This is a workaround until MDAnalysis implements its own CIF reader.
"""

import MDAnalysis as mda
from openbabel import openbabel as ob


def convert_cif_to_pdb(cif_file: str):
    """Convert a CIF file to a PDB file. Uses OpenBabel to do the conversion."""
    ob_conv = ob.OBConversion()
    ob_conv.SetInFormat("cif")
    mol = ob.OBMol()
    ob_conv.ReadFile(mol, cif_file)
    ob_conv.SetOutFormat("pdb")
    pdb_str = ob_conv.WriteString(mol)
    return pdb_str, mol
