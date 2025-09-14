"""
Exposes the compiled RDKit bindings under `lahuta.rdkit`, so imports
like `from lahuta.rdkit import RWMol` work and static analyzers can
resolve the module.

At runtime, the module simply re-exports everything from `lahuta.lib.lahuta.rdkit`
"""

from __future__ import annotations

# fmt: off
try:
    # Import the main lahuta module and re-export its rdkit submodule
    from .lib import lahuta as lxx

    # Re-export all rdkit symbols explicitly
    _rdkit = lxx.rdkit

    # Export all rdkit classes and functions
    Atom       = _rdkit.Atom
    Bond       = _rdkit.Bond
    RWMol      = _rdkit.RWMol
    Point3D    = _rdkit.Point3D
    BondDir    = _rdkit.BondDir
    BondStereo = _rdkit.BondStereo
    BondType   = _rdkit.BondType
    Conformer  = _rdkit.Conformer
    AtomMonomerInfo     = _rdkit.AtomMonomerInfo
    AtomPDBResidueInfo  = _rdkit.AtomPDBResidueInfo
    pAtomPDBResidueInfo = _rdkit.pAtomPDBResidueInfo
    hasNonZeroZCoords   = _rdkit.hasNonZeroZCoords

    AROMATIC       = _rdkit.AROMATIC
    BEGINDASH      = _rdkit.BEGINDASH
    BEGINWEDGE     = _rdkit.BEGINWEDGE
    DOUBLE         = _rdkit.DOUBLE
    EITHERDOUBLE   = _rdkit.EITHERDOUBLE
    ENDDOWNRIGHT   = _rdkit.ENDDOWNRIGHT
    ENDUPRIGHT     = _rdkit.ENDUPRIGHT
    NONE           = _rdkit.NONE
    SINGLE         = _rdkit.SINGLE
    STEREOANY      = _rdkit.STEREOANY
    STEREOATROPCCW = _rdkit.STEREOATROPCCW
    STEREOATROPCW  = _rdkit.STEREOATROPCW
    STEREOCIS      = _rdkit.STEREOCIS
    STEREOE        = _rdkit.STEREOE
    STEREONONE     = _rdkit.STEREONONE
    STEREOTRANS    = _rdkit.STEREOTRANS
    STEREOZ        = _rdkit.STEREOZ
    TRIPLE         = _rdkit.TRIPLE
    UNKNOWN        = _rdkit.UNKNOWN
    UNSPECIFIED    = _rdkit.UNSPECIFIED

except Exception as _e:
    raise ImportError(
        "lahuta.rdkit requires the compiled C++ bindings. "
        "Build the core and reinstall the Python package. "
        f"Original error: {_e}"
    ) from _e
