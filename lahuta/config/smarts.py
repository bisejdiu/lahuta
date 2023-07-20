from enum import Enum

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


_Atom_Type_Names: list[str] = [member.name for member in SmartsPatternRegistry]
available_atom_types_dict = {name: i for i, name in enumerate(_Atom_Type_Names)}
AVAILABLE_ATOM_TYPES = Enum("AVAILABLE_ATOM_TYPES", available_atom_types_dict)
