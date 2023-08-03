"""
Module: smarts.py

This module defines a registry for SMARTS patterns associated with different atom types, and maps 
atom types to unique integer identifiers.

Classes:
    SmartsPatternRegistry(Enum): Enum class that serves as a registry for SMARTS patterns 
                                 associated with different atom types. Each atom type is represented
                                 by a string key and a corresponding SMARTS pattern.

Variables:
    _Atom_Type_Names (list): A list of names of members in the SmartsPatternRegistry.

    AVAILABLE_ATOM_TYPES (dict): Dictionary that maps each atom type to a unique integer identifier.
                                 The atom types are strings and the identifiers are indices from 
                                 _Atom_Type_Names.

Notes:
    SMARTS (SMiles ARbitrary Target Specification) is a language that allows you to specify 
    substructural patterns in molecules. It is an extension of the simpler SMILES notation. 
    Each SMARTS pattern in the registry represents a specific type of atom or chemical group.
"""

from enum import Enum
from typing import Dict, List

import lahuta.config._smart_strings as smarts


# "NH2 terminal amide" and "oxygen amide term" are ambiguous and may affect non-protein amides
class SmartsPatternRegistry(Enum):
    """A registry of SMARTS patterns for atom typing.

    Each atom type is represented by a string key and a SMARTS pattern.
    """

    HBOND_ACCEPTOR = smarts.HBOND_ACCEPTOR_SMARTS
    HBOND_DONOR = smarts.HBOND_DONOR_SMARTS

    # nonsensical string to avoid duplicate entries in Enum
    XBOND_ACCEPTOR = {**smarts.HBOND_ACCEPTOR_SMARTS, "": "[Xx]"}

    XBOND_DONOR = {"donor": smarts.SMARTS_STR_XBD}

    WEAK_HBOND_ACCEPTOR = {
        **smarts.HBOND_ACCEPTOR_SMARTS,
        "c-x halogens": smarts.SMARTS_STR_XBD,
    }

    WEAK_HBOND_DONOR = {"donor": smarts.SMARTS_STR_WHBD}
    POS_IONISABLE = smarts.POSITIVELY_IOINISABLE_SMARTS

    NEG_IONISABLE = {
        "O acidic group": smarts.SMARTS_STR_NI_O_ACIDIC_GRP,
        "anions": smarts.SMARTS_STR_NI_ANIONS,
    }

    HYDROPHOBE = {"hydrophobe": smarts.SMARTS_STR_HYDROPHOBIC}
    CARBONYL_OXYGEN = {"oxygen": smarts.SMARTS_STR_CARBONYL_OXYGEN}
    CARBONYL_CARBON = {"carbon": smarts.SMARTS_STR_CARBONYL_CARBON}
    AROMATIC = smarts.AROMATIC_SMARTS


_Atom_Type_Names: List[str] = [member.name for member in SmartsPatternRegistry]
AVAILABLE_ATOM_TYPES: Dict[str, int] = {name: i for i, name in enumerate(_Atom_Type_Names)}
