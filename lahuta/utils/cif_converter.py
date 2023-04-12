"""
Converts CIF files to PDB files. This is a workaround until MDAnalysis implements its own CIF reader.
"""

import tempfile

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


def read_cif_as_pdb(cif_file: str):
    """Convert a CIF file to a PDB file. Uses MDAnalysis to do the conversion."""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as pdb_file:

        ob_conv = ob.OBConversion()
        ob_conv.SetInFormat("cif")
        temp_mol = ob.OBMol()
        ob_conv.ReadFile(temp_mol, cif_file)

        ob_conv.SetOutFormat("pdb")
        ob_conv.WriteFile(temp_mol, pdb_file.name)
        pdb_str = ob_conv.WriteString(temp_mol)

        pdb_file.write(pdb_str)
        pdb_file.flush()

        ob_conv = ob.OBConversion()
        ob_conv.SetInFormat("pdb")
        mol = ob.OBMol()
        ob_conv.ReadString(mol, pdb_str)

        return pdb_str, mol
