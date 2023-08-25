"""Provides definitions for various sets of SMARTS patterns used in to find and match
specific atom types. Additionally, it includes a helper function for remapping atom types for
easier lookups.

"""

from enum import Enum
from typing import ClassVar

import lahuta.config._smart_strings as smarts


# "NH2 terminal amide" and "oxygen amide term" are ambiguous and may affect non-protein amides
class SmartsPatternRegistry(Enum):
    """A registry of SMARTS patterns for atom typing.

    Each atom type is represented by a string key and a SMARTS pattern. See the
    definitions below for further details, and follow the links for the corresponding
    SMARTS patterns.
    """

    HBOND_ACCEPTOR = smarts.HBOND_ACCEPTOR_SMARTS
    """
    Hydrogen bond acceptor SMARTS patterns. See 
    [definitions](smarts_defs.md#lahuta.config._smart_strings.HBOND_ACCEPTOR_SMARTS) for details.
    """
    HBOND_DONOR = smarts.HBOND_DONOR_SMARTS

    # nonsensical string to avoid duplicate entries in Enum
    XBOND_ACCEPTOR: ClassVar = {**smarts.HBOND_ACCEPTOR_SMARTS, "": "[Xx]"}
    """
    Halogen bond acceptor SMARTS patterns. See
    [definitions](smarts_defs.md#lahuta.config._smart_strings.HBOND_ACCEPTOR_SMARTS) for details.
    Note that we provide a nonsensical string to avoid duplicate entries in the Enum.
    """

    XBOND_DONOR: ClassVar = {"donor": smarts.SMARTS_STR_XBD}
    """
    Halogen bond donor SMARTS patterns. See
    [definitions](smarts_defs.md#lahuta.config._smart_strings.SMARTS_STR_XBD) for details.
    """

    WEAK_HBOND_ACCEPTOR: ClassVar = {
        **smarts.HBOND_ACCEPTOR_SMARTS,
        "c-x halogens": smarts.SMARTS_STR_XBD,
    }
    """
    Weak hydrogen bond acceptor SMARTS patterns. In addition to the standard
    [hydrogen bond acceptor](smarts_defs.md#lahuta.config._smart_strings.HBOND_ACCEPTOR_SMARTS)
    SMARTS patterns, this also includes
    [halogen bond acceptor](smarts_defs.md#lahuta.config._smart_strings.HBOND_ACCEPTOR_SMARTS) 
    SMARTS patterns.
    """

    WEAK_HBOND_DONOR: ClassVar = {"donor": smarts.SMARTS_STR_WHBD}
    """
    Weak hydrogen bond donor SMARTS patterns. See
    [definitions](smarts_defs.md#lahuta.config._smart_strings.SMARTS_STR_WHBD) for details.
    """

    POS_IONISABLE: ClassVar = smarts.POSITIVELY_IOINISABLE_SMARTS
    """
    Positively ionisable SMARTS patterns. See
    [definitions](smarts_defs.md#lahuta.config._smart_strings.POSITIVELY_IOINISABLE_SMARTS) for
    details.
    """

    NEG_IONISABLE: ClassVar = {
        "O acidic group": smarts.SMARTS_STR_NI_O_ACIDIC_GRP,
        "anions": smarts.SMARTS_STR_NI_ANIONS,
    }
    """
    Negatively ionisable SMARTS patterns. It includes two types of patterns: 
    [acidic oxygen](smarts_defs.md#lahuta.config._smart_strings.SMARTS_STR_NI_O_ACIDIC_GRP) and
    [anions](smarts_defs.md#lahuta.config._smart_strings.SMARTS_STR_NI_ANIONS).
    """

    HYDROPHOBE: ClassVar = {"hydrophobe": smarts.SMARTS_STR_HYDROPHOBIC}
    """
    Hydrophobic SMARTS patterns. See
    [definitions](smarts_defs.md#lahuta.config._smart_strings.SMARTS_STR_HYDROPHOBIC) for details.
    """

    CARBONYL_OXYGEN: ClassVar = {"oxygen": smarts.SMARTS_STR_CARBONYL_OXYGEN}
    """
    Carbonyl oxygen SMARTS patterns. See
    [definitions](smarts_defs.md#lahuta.config._smart_strings.SMARTS_STR_CARBONYL_OXYGEN) for
    details.
    """

    CARBONYL_CARBON: ClassVar = {"carbon": smarts.SMARTS_STR_CARBONYL_CARBON}
    """
    Carbonyl carbon SMARTS patterns. See
    [definitions](smarts_defs.md#lahuta.config._smart_strings.SMARTS_STR_CARBONYL_CARBON) for
    details.
    """

    AROMATIC = smarts.AROMATIC_SMARTS
    """
    Aromatic SMARTS patterns. See
    [definitions](smarts_defs.md#lahuta.config._smart_strings.AROMATIC_SMARTS) for details.
    """


_Atom_Type_Names: list[str] = [member.name for member in SmartsPatternRegistry]
AVAILABLE_ATOM_TYPES: dict[str, int] = {name: i for i, name in enumerate(_Atom_Type_Names)}
"""
Type: `dict[str, int]`: A dictionary mapping atom type names (strings) to unique integer identifiers.
"""
