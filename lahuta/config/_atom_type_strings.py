"""
module: `config._atom_type_strings.py`

This module provides string representations of various atom types and categories.

The various types of atoms represented here include metals, standard amino acids, 
different types of bond acceptors and donors, ionisable atoms, 
hydrophobic atoms, carbonyl atoms, and aromatic atoms.

Variables:
    ```
    _METALS_STR (str): A string of metal atom types.
    _STANDARD_AA_STR (str): A string of standard amino acid atom types.
    _HA_ATOM_TYPES (str): A string of hydrogen bond acceptor atom types.
    _HD_ATOM_TYPES (str): A string of hydrogen bond donor atom types.
    _XA_ATOM_TYPES (str): A string of halogen bond acceptor atom types.
    _XD_ATOM_TYPES (str): A string of halogen bond donor atom types.
    _WHA_ATOM_TYPES (str): A string of weak hydrogen bond acceptor atom types.
    _WHD_ATOM_TYPES (str): A string of weak hydrogen bond donor atom types.
    _POS_IONISABLE_ATOM_TYPES (str): A string of positive ionisable atom types.
    _NEG_IONISABLE_ATOM_TYPES (str): A string of negative ionisable atom types.
    _HYDROPHOBE_ATOM_TYPES (str): A string of hydrophobic atom types.
    _CARBONYL_OXYGEN_ATOM_TYPES (str): A string of carbonyl oxygen atom types.
    _CARBONYL_CARBON_ATOM_TYPES (str): A string of carbonyl carbon atom types.
    _AROMATIC_ATOM_TYPES (str): A string of aromatic atom types.
    RESIDUE_SYNONYMS (dict): A dictionary of residue names and their synonyms.
    STANDARD_AMINO_ACIDS (set): A set of standard amino acids.

    ```

The atom types are initially defined as comma-separated strings. 
These strings are then processed into sets for easy and efficient access 
throughout the rest of the library.

"""
from typing import Dict, List, Set

from MDAnalysis.core.selection import ProteinSelection

# fmt: off
_METALS_STR = (
    "Li,Be,Na,Mg,Aa,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Rb,Sr,Y,Zr,Nb,Mo,"
    "Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,"
    "Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf"
)

# _STANDARD_AA_STR = "ALA,CYS,ASP,GLU,PHE,GLY,HIS,HSE,HSD,ILE,LYS,LEU,MET,ASN,PRO,GLN,ARG,SER,THR,VAL,TRP,TYR"

_HA_ATOM_TYPES = (
    "ALAO,ARGO,ASNO,ASPO,CYSO,GLNO,GLUO,GLYO,HISO,ILEO,LEUO,LYSO,METO,PHEO,PROO,"
    "SERO,THRO,TRPO,TYRO,VALO,ALAOXT,ARGOXT,ASNOXT,ASPOXT,CYSOXT,GLNOXT,GLUOXT,GLYOXT,"
    "HISOXT,ILEOXT,LEUOXT,LYSOXT,METOXT,PHEOXT,PROOXT,SEROXT,THROXT,TRPOXT,TYROXT,VALOXT,"
    "ASNOD1,ASNND2,ASPOD1,ASPOD2,GLNOE1,GLNNE2,GLUOE1,GLUOE2,HISND1,HISCE1,HISNE2,"
    "HISCD2,METSD,CYSSG,SEROG,THROG1,TYROH"
)

_HD_ATOM_TYPES = (
    "ALAN,ARGN,ASNN,ASPN,CYSN,GLNN,GLUN,GLYN,HISN,ILEN,LEUN,LYSN,METN,PHEN,"
    "SERN,THRN,TRPN,TYRN,VALN,ARGNE,ARGNH1,ARGNH2,ASNND2,ASNOD1,CYSSG,GLNNE2,"
    "GLNOE1,HISND1,HISCE1,HISNE2,HISCD2,LYSNZ,SEROG,THROG1,TRPNE1,TYROH"
)

_XA_ATOM_TYPES = (
    "ALAO,ARGO,ASNO,ASPO,CYSO,GLNO,GLUO,GLYO,HISO,ILEO,LEUO,LYSO,METO,PHEO,PROO,"
    "SERO,THRO,TRPO,TYRO,VALO,ALAOXT,ARGOXT,ASNOXT,ASPOXT,CYSOXT,GLNOXT,GLUOXT,"
    "GLYOXT,HISOXT,ILEOXT,LEUOXT,LYSOXT,METOXT,PHEOXT,PROOXT,SEROXT,THROXT,TRPOXT,"
    "TYROXT,VALOXT,ASNOD1,ASNND2,ASPOD1,ASPOD2,GLNOE1,GLNNE2,GLUOE1,GLUOE2,HISND1,"
    "HISCE1,HISNE2,HISCD2,METSD,CYSSG,SEROG,THROG1,TYROH"
)

_XD_ATOM_TYPES = ""

_WHA_ATOM_TYPES = (
    "ALAO,ARGO,ASNO,ASPO,CYSO,GLNO,GLUO,GLYO,HISO,ILEO,LEUO,LYSO,METO,PHEO,"
    "PROO,SERO,THRO,TRPO,TYRO,VALO,ALAOXT,ARGOXT,ASNOXT,ASPOXT,CYSOXT,GLNOXT,"
    "GLUOXT,GLYOXT,HISOXT,ILEOXT,LEUOXT,LYSOXT,METOXT,PHEOXT,PROOXT,SEROXT,THROXT,"
    "TRPOXT,TYROXT,VALOXT,ASNOD1,ASNND2,ASPOD1,ASPOD2,GLNOE1,GLNNE2,GLUOE1,GLUOE2,"
    "HISND1,HISCE1,HISNE2,HISCD2,METSD,CYSSG,SEROG,THROG1,TYROH"
)

_WHD_ATOM_TYPES = (
    "ALACA,ARGCA,ASNCA,ASPCA,CYSCA,GLNCA,GLUCA,GLYCA,HISCA,ILECA,LEUCA,LYSCA,"
    "METCA,PHECA,PROCA,SERCA,THRCA,TRPCA,TYRCA,VALCA,ALACB,ARGCB,ARGCG,ARGCD,"
    "ASNCB,ASPCB,CYSCB,GLNCB,GLNCG,GLUCB,GLUCG,GLNCB,HISCB,ILECB,ILECG1,ILECD1,"
    "ILECG2,LEUCB,LEUCG,LEUCD1,LEUCD2,LYSCB,LYSCG,LYSCD,LYSCE,METCB,METCG,METCE,"
    "PHECB,PHECG,PHECD1,PHECD2,PHECE1,PHECE2,PHECZ,PROCB,PROCG,PROCD,SERCB,"
    "THRCB,THRCG2,TRPCB,TRPCD1TRPCE3,TRPCZ3,TRPCH2,TRPCZ2,TYRCB,TYRCD1,TYRCD2,"
    "TYRCE1,TYRCE2,TRYCB,VALCB,VALCG1,VALCG2"
)

_POS_IONISABLE_ATOM_TYPES = "ARGNE,ARGCZ,ARGNH1,ARGNH2,HISCG,HISND1,HISCE1,HISNE2,HISCD2,LYSNZ"

_NEG_IONISABLE_ATOM_TYPES = "ASPOD1,ASPOD2,GLUOE1,GLUOE2"

_HYDROPHOBE_ATOM_TYPES = (
    "ALACB,ARGCB,ARGCG,ASNCB,ASPCB,CYSCB,GLNCB,GLNCG,GLUCB,GLUCG,GLNCB,HISCB,"
    "ILECB,ILECG1,ILECD1,ILECG2,LEUCB,LEUCG,LEUCD1,LEUCD2,LYSCB,LYSCG,LYSCD,"
    "METCB,METCG,METSD,METCE,PHECB,PHECG,PHECD1,PHECD2,PHECE1,PHECE2,PHECZ,PROCB,"
    "PROCG,THRCG2,TRPCB,TRPCG,TRPCD2,TRPCE3,TRPCZ3,TRPCH2,TRPCZ2,TRYCB,TYRCG,"
    "TYRCD1,TYRCD2,TYRCE1,TYRCE2,VALCB,VALCG1,VALCG2"
)

_CARBONYL_OXYGEN_ATOM_TYPES = (
    "ALAO,ARGO,ASNO,ASPO,CYSO,GLNO,GLUO,GLYO,HISO,ILEO,LEUO," 
    "LYSO,METO,PHEO,PROO,SERO,THRO,TRPO,TYRO,VALO"
)

_CARBONYL_CARBON_ATOM_TYPES = (
    "ALAC,ARGC,ASNC,ASPC,CYSC,GLNC,GLUC,GLYC,HISC,ILEC,LEUC,LYSC," 
    "METC,PHEC,PROC,SERC,THRC,TRPC,TYRC,VALC"
)

_AROMATIC_ATOM_TYPES = (
    "HISCG,HISND1,HISCE1,HISNE2,HISCD2,PHECG,PHECD1,PHECD2,PHECE1,PHECE2,PHECZ,"
    "TRPCG,TRPCD1,TRPCD2,TRPNE1,TRPCE2,TRPCE3,TRPCZ2,TRPCZ3,TRPCH2,TYRCG,TYRCD1,"
    "TYRCD2,TYRCE1,TYRCE2,TYRCZ"
)

RESIDUE_SYNONYMS = {
    "HIS": [
        "HIS", "HSD", "HSE", "HSP", "HIE", "HIP", "HID", "HIS1", 
        "HIS2", "HISA", "HISB", "HISD", "HISE", "HISH", "HYP"
    ],
    "PHE": ["PHE"],
    "TYR": ["TYR"],
    "ALA": ["ALA", "ALAD"],
    "ARG": ["ARG", "ARGN"],
    "ASN": ["ASN", "ASN1"],
    "ASP": ["ASP", "ASPH"],
    "CYS": ["CYS", "CYS1", "CYS2", "CYSH", "CYX", "CYM"],
    "GLN": ["GLN", "GLH"],
    "GLU": ["GLU", "GLUH"],
    "GLY": ["GLY"],
    "ILE": ["ILE"],
    "LEU": ["LEU"],
    "LYS": ["LYS", "LYSH", "LYN"],
    "MET": ["MET"],
    "PRO": ["PRO"],
    "SER": ["SER"],
    "THR": ["THR"],
    "TRP": ["TRP"],
    "VAL": ["VAL"],
}

BASE_AA_CONVERSION: Dict[str, str] = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLN': 'Q',
    'GLU': 'E',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V'
}

def parse_atom_types_string(_atom_types_string: str) -> Dict[str, List[str]]:
    """
    Parse a string of atom types into a dictionary of residue names and atom parts.

    Args:
        _atom_types_string (str): A string of atom types.

    Returns:
        Dict[str, List[str]]: A dictionary of residue names and atom parts.
    """

    atom_parts: Dict[str, List[str]] = {}
    for atom_type in _atom_types_string.split(","):
        residue_name = atom_type[:3]
        atom_part = atom_type[3:]

        if residue_name not in atom_parts:
            atom_parts[residue_name] = []
        atom_parts[residue_name].append(atom_part)

    return atom_parts

def parse_atom_types(_atom_types_string: str) -> Set[str]:
    """
    Parse a string of atom types into a set of atom types.

    Args:
        _atom_types_string (str): A string of atom types.

    Returns:
        Set[str]: A set of atom types.
    """

    res_atoms = parse_atom_types_string(_atom_types_string)

    atom_types: Set[str] = set()
    for residue, synonyms in RESIDUE_SYNONYMS.items():
        for synonym in synonyms:
            atom_parts = res_atoms.get(residue)
            if atom_parts is None:
                continue
            for atom_part in atom_parts:
                atom_types.add(synonym + atom_part)

    return atom_types

# FIXME: update with the residue list provided by MDAnalysis
METALS = set(_METALS_STR.split(","))
# STANDARD_AMINO_ACIDS = set(_STANDARD_AA_STR.split(","))
STANDARD_AMINO_ACIDS = ProteinSelection.prot_res
"""
module: `lahuta.config._atom_type_strings.py`

Type: `Set[str]`: A set of standard amino acids. Taken from `MDAnalysis.core.selection.ProteinSelection`.
"""
HBOND_ACCEPTORS = parse_atom_types(_HA_ATOM_TYPES) # set(_HA_ATOM_TYPES.split(","))
HBOND_DONORS = parse_atom_types(_HD_ATOM_TYPES) #set(_HD_ATOM_TYPES.split(","))
XBOND_ACCEPTORS = parse_atom_types(_XA_ATOM_TYPES) #set(_XA_ATOM_TYPES.split(","))
XBOND_DONORS = set("") #set("")
WEAK_HBOND_ACCEPTORS = parse_atom_types(_WHA_ATOM_TYPES) #set(_WHA_ATOM_TYPES.split(","))
WEAK_HBOND_DONORS = parse_atom_types(_WHD_ATOM_TYPES) #set(_WHD_ATOM_TYPES.split(","))
POS_IONISABLE = parse_atom_types(_POS_IONISABLE_ATOM_TYPES) #set(_POS_IONISABLE_ATOM_TYPES.split(","))
NEG_IONISABLE = parse_atom_types(_NEG_IONISABLE_ATOM_TYPES) #set(_NEG_IONISABLE_ATOM_TYPES.split(","))
HYDROPHOBES = parse_atom_types(_HYDROPHOBE_ATOM_TYPES) #set(_HYDROPHOBE_ATOM_TYPES.split(","))
CARBONYL_OXYGENS = parse_atom_types(_CARBONYL_OXYGEN_ATOM_TYPES) #set(_CARBONYL_OXYGEN_ATOM_TYPES.split(","))
CARBONYL_CARBONS = parse_atom_types(_CARBONYL_CARBON_ATOM_TYPES) #set(_CARBONYL_CARBON_ATOM_TYPES.split(","))
AROMATIC = parse_atom_types(_AROMATIC_ATOM_TYPES) #set(_AROMATIC_ATOM_TYPES.split(","))

# NOTES
# "ALAO",    # all the carbonyl Oxygens in the main chain
# "ALAOXT",  # all the carbonyl Oxygens terminals
# "ASNND2",  # for the ambiguity of the position of the N and O
# "GLNOE1"   # for the ambiguity of the position of the N and O
# "GLNNE2",  # for the ambiguity of the position of the N and O
# "HISND1",  # for the ambiguity of the position of the N/C
# "HISCE1",  # for the ambiguity of the position of the N/C
# "HISNE2",  # for the ambiguity of the position of the N/C
# "HISCD2",  # for the ambiguity of the position of the N/C
# "METSD",   # http://pubs.acs.org/doi/abs/10.1021/jz300207k and pubid 19089987
# "CYSSG",   # pubid 19089987, also when they from di-sulfide (Cys-Cys, fig 8 paper)
# "SEROG",   # isostar plots
# "THROG1",  # isostar plots
# "TYROH",   # isostar plots
# "ALAN",    # all the amide nitrogens in the main chain except proline
# "ASNOD1",  # for the ambiguity of the position of the N and O
# "CYSSG",   # http://www.ncbi.nlm.nih.gov/pubmed/19089987
# "ALACA",   # all the c-alphas
# "ALACB",   # cb and further down
# "CYSCB",   # sulfur in Cys has an Hydrogen, it is polarised
