"""Contains a collection of Simplified Molecular-Input Line-Entry System
(SMARTS) strings used for SMARTS pattern matching.

The SMARTS strings cater to a variety of chemical functionalities and patterns,
including but not limited to hydrogen bond acceptors and donors, positively and negatively
ionizable groups, hydrophobic functionalities, carbonyl groups, and aromatic systems.

Each SMARTS string represents a specific chemical pattern or functionality, with separate 
SMARTS strings defined for different types of hydrogen bond acceptors, such as a standard 
acceptor, enol, tautomeric NH, and terminal amide.

The SMARTS strings are grouped into dictionaries based on their functionality. For example,
HBOND_ACCEPTOR_SMARTS includes all SMARTS strings for different types of hydrogen bond acceptors,
while POSITIVELY_IOINISABLE_SMARTS contains SMARTS strings representing different types of
positively ionizable groups.

The module also defines aromatic systems of varying ring sizes, ranging from 4 to 8.

Notes:
    SMARTS is a language used to specify substructural patterns in molecules. It extends the
    simpler SMILES syntax and provides a more flexible way to describe molecular patterns.

    This module does not contain any functions or classes, but only data in the form of strings and
    dictionaries of strings.

    It is expected to be used as a source of SMARTS strings for other modules and applications
    dealing with molecular structure processing, matching, and analysis.

"""
# HYDROGEN BOND ACCEPTOR
SMARTS_STR_HBA_ACCEPTOR = (
    "[#8,#9,$([#16;H0,H1;v2,v1]),$([N;v3;!$(N-*=!@[O,N,P,S]);" "!$(N-!@a);!$([NH]=!@*)]),$([nH0;+0])]"
)
"""
Type: `str`: Standard hydrogen bond acceptor. Helps define 
[HBOND_ACCEPTOR_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.HBOND_ACCEPTOR_SMARTS).
"""
SMARTS_STR_HBA_ENOL = "[$([nH]:@c(=O))]"
"""
Type: `str`: Enol. Helps define 
[HBOND_ACCEPTOR_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.HBOND_ACCEPTOR_SMARTS).
"""
SMARTS_STR_HBA_TAUTOMERIC_NH = "[$([n;H1;v3;!$([nH]cccc)])]"
"""
Type: `str`: Tautomeric NH. Helps define
[HBOND_ACCEPTOR_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.HBOND_ACCEPTOR_SMARTS).
"""
SMARTS_STR_HBA_TERMINAL_AMIDE = "[$([N;H2;v3;$(N-C(=O))])]"
"""
Type: `str`: NH2 terminal amide. Helps define
[HBOND_ACCEPTOR_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.HBOND_ACCEPTOR_SMARTS).
"""
# HYDROGEN BOND DONOR
SMARTS_STR_HBD_DONOR = "[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]"
"""
Type: `str`: Standard hydrogen bond donor. Helps define
[HBOND_DONOR_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.HBOND_DONOR_SMARTS).
"""
SMARTS_STR_HBD_OXYGEN_ACID = "[$([O;H0;$(O=C([OH])-*)])]"
"""
Type: `str`: Oxygen acid. Helps define
[HBOND_DONOR_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.HBOND_DONOR_SMARTS).
"""
SMARTS_STR_HBD_TAUTOMERIC_NH = "[$(n:a:[nH])]"
"""
Type: `str`: Tautomeric NH. Helps define
[HBOND_DONOR_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.HBOND_DONOR_SMARTS).
"""
SMARTS_STR_HBD_OXYGEN_AMIDE_TERM = "[$([O;H0;$(O=C-[NH2])])]"
"""
Type: `str`: Oxygen amide terminal. Helps define
[HBOND_DONOR_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.HBOND_DONOR_SMARTS).
"""
# WEAK HBOND DONOR and HALOGEN BOND DONOR
SMARTS_STR_WHBD = "[#6!H0]"
"""
Type: `str`: Weak hydrogen bond donor. Helps define
[WEAK_HBOND_DONOR](smarts_defs.md#lahuta.config.smarts_patterns.WEAK_HBOND_DONOR).
"""

SMARTS_STR_XBD = "[Cl,Br,I;X1;$([Cl,Br,I]-[#6])]"
"""
Type: `str`: Halogen bond donor. Helps define
[XBOND_DONOR](smarts_defs.md#lahuta.config.smarts_patterns.XBOND_DONOR).
"""

# POSITIVELY IONISABLE
SMARTS_STR_PI_RDKIT_BASIC_GRP = (
    "[$([N;H2&+0][C;!$(C=*)]),$([N;H1&+0]([C;!$(C=*)])[C;!$(C=*)]),"
    "$([N;H0&+0]([C;!$(C=*)])([C;!$(C=*)])[C;!$(C=*)]);!$(N[a])]"
)
"""
Type: `str`: RDKIT basic group. Helps define
[POSITIVELY_IOINISABLE_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.POSITIVELY_IOINISABLE_SMARTS).
"""
SMARTS_STR_PI_IMIDAZOLE = "[n;R1]1[c;R1][n;R1][c;R1][c;R1]1"
"""
Type: `str`: Imidazole. Helps define
[POSITIVELY_IOINISABLE_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.POSITIVELY_IOINISABLE_SMARTS).
"""
SMARTS_STR_PI_GUANIDINE_AMIDINE = "NC(=N)"
"""
Type: `str`: Guanidine amidine. Helps define
[POSITIVELY_IOINISABLE_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.POSITIVELY_IOINISABLE_SMARTS).
"""
SMARTS_STR_PI_RDKIT_POSN = "[#7;+;!$([N+]-[O-])]"
"""
Type: `str`: RDKIT positive nitrogen. Helps define
[POSITIVELY_IOINISABLE_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.POSITIVELY_IOINISABLE_SMARTS).
"""
SMARTS_STR_PI_CATIONS = "[$([*+1,*+2,*+3]);!$([N+]-[O-])]"
"""
Type: `str`: Cations. Helps define
[POSITIVELY_IOINISABLE_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.POSITIVELY_IOINISABLE_SMARTS).
"""
SMARTS_STR_PI_METALS = (
    "[Li,Be,Na,Mg,Al,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Rb,Sr,Y,Zr,"
    "Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,"
    "Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,Fr,Ra,Ac,Th,"
    "Pa,U,Np,Pu,Am,Cm,Bk,Cf]"
)
"""
Type: `str`: Metals. Helps define
[POSITIVELY_IOINISABLE_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.POSITIVELY_IOINISABLE_SMARTS).
"""
# NEGATIVELY IONISABLE
SMARTS_STR_NI_ANIONS = "[*-1,*-2]"
"""
Type: `str`: Anions. Helps define
[NEG_IONISABLE](smarts_defs.md#lahuta.config.smarts_patterns.NEG_IONISABLE).
"""
SMARTS_STR_NI_O_ACIDIC_GRP = "[$([OH,O-]-[C,S,N,P,Cl,Br,I]=O),$(O=[C,S,N,P,Cl,Br,I]-[OH,O-])]"
"""
Type: `str`: Oxygen acidic group. Helps define
[NEG_IONISABLE](smarts_defs.md#lahuta.config.smarts_patterns.NEG_IONISABLE).
"""

# HYDROPHOBIC, and CARBONYL
SMARTS_STR_HYDROPHOBIC = "[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,Cl+0,Br+0,I+0]"
"""
Type: `str`: Hydrophobic. Helps define
[HYDROPHOBE](smarts_defs.md#lahuta.config.smarts_patterns.HYDROPHOBE).
"""
SMARTS_STR_CARBONYL_OXYGEN = "[$([OH0]=[CX3,c]);!$([OH0]=[CX3,c]-[OH,O-])]"
"""
Type: `str`: Carbonyl oxygen. Helps define
[CARBONYL_OXYGEN](smarts_defs.md#lahuta.config.smarts_patterns.CARBONYL_OXYGEN).
"""
SMARTS_STR_CARBONYL_CARBON = "[$([CX3,c]=[OH0]);!$([CX3,c](=[OH0])-[OH,O-])]"
"""
Type: `str`: Carbonyl carbon. Helps define
[CARBONYL_CARBON](smarts_defs.md#lahuta.config.smarts_patterns.CARBONYL_CARBON).
"""

# AROMATIC
SMARTS_STR_AROMATIC_4 = "[a;r4,!R1&r3]1:[a;r4,!R1&r3]:[a;r4,!R1&r3]:[a;r4,!R1&r3]:1"
"""
Type: `str`: Aromatic 4. Helps define
[AROMATIC_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.AROMATIC_SMARTS).
"""
SMARTS_STR_AROMATIC_5 = (
    "[a;r5,!R1&r4,!R1&r3]1:[a;r5,!R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:" "[a;r5,!R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:1"
)
"""
Type: `str`: Aromatic 5. Helps define
[AROMATIC_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.AROMATIC_SMARTS).
"""
SMARTS_STR_AROMATIC_6 = (
    "[a;r6,!R1&r5,!R1&r4,!R1&r3]1:[a;r6,!R1&r5,!R1&r4,!R1&r3]:"
    "[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:"
    "[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:1"
)
"""
Type: `str`: Aromatic 6. Helps define
[AROMATIC_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.AROMATIC_SMARTS).
"""
SMARTS_STR_AROMATIC_7 = (
    "[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]1:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:"
    "[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:"
    "[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:"
    "[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:1"
)
"""
Type: `str`: Aromatic 7. Helps define
[AROMATIC_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.AROMATIC_SMARTS).
"""
SMARTS_STR_AROMATIC_8 = (
    "[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]1:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:"
    "[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:"
    "[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:"
    "[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:1"
)
"""
Type: `str`: Aromatic 8. Helps define
[AROMATIC_SMARTS](smarts_defs.md#lahuta.config.smarts_patterns.AROMATIC_SMARTS).
"""

# AMIDE, DEFINITION FROM `http://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html`
AMIDE_SMARTS = "[NX3][CX3](=[OX1])[#6]"
HBOND_ACCEPTOR_SMARTS = {
    "acceptor": SMARTS_STR_HBA_ACCEPTOR,
    "enol": SMARTS_STR_HBA_ENOL,
    "tautomeric nH": SMARTS_STR_HBA_TAUTOMERIC_NH,
    "NH2 terminal amide": SMARTS_STR_HBA_TERMINAL_AMIDE,
}
"""
A dictionary of hydrogen bond acceptor SMARTS patterns.
!!! tip "Definitions"
    - `acceptor`: Standard hydrogen bond acceptor \
        ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_HBA_ACCEPTOR))
    - `enol`: Enol ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_HBA_ENOL))
    - `tautomeric nH`: Tautomeric NH \
        ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_HBA_TAUTOMERIC_NH))
    - `NH2 terminal amide`: NH2 terminal amide \
        ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_HBA_TERMINAL_AMIDE))
"""

POSITIVELY_IOINISABLE_SMARTS = {
    "rdkit basic group": SMARTS_STR_PI_RDKIT_BASIC_GRP,
    "imidazole": SMARTS_STR_PI_IMIDAZOLE,
    "guanidine amidine": SMARTS_STR_PI_GUANIDINE_AMIDINE,
    "rdkit posn": SMARTS_STR_PI_RDKIT_POSN,
    "cations": SMARTS_STR_PI_CATIONS,
    "metals": SMARTS_STR_PI_METALS,
}
"""
A dictionary of positively ionizable SMARTS patterns.
!!! tip "Definitions"
    - `rdkit basic group`: RDKIT basic group \
        ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_PI_RDKIT_BASIC_GRP))
    - `imidazole`: Imidazole ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_PI_IMIDAZOLE))
    - `guanidine amidine`: Guanidine amidine \
        ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_PI_GUANIDINE_AMIDINE))
    - `rdkit posn`: RDKIT positive nitrogen \
        ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_PI_RDKIT_POSN))
    - `cations`: Cations ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_PI_CATIONS))
    - `metals`: Metals ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_PI_METALS))
"""

HBOND_DONOR_SMARTS = {
    "donor": SMARTS_STR_HBD_DONOR,
    "oxygen acid": SMARTS_STR_HBD_OXYGEN_ACID,
    "tautomer nH": SMARTS_STR_HBD_TAUTOMERIC_NH,
    "oxygen amide term": SMARTS_STR_HBD_OXYGEN_AMIDE_TERM,
}
"""
A dictionary of hydrogen bond donor SMARTS patterns.
!!! tip "Definitions"
    - `donor`: Standard hydrogen bond donor ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_HBD_DONOR))
    - `oxygen acid`: Oxygen acid ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_HBD_OXYGEN_ACID))
    - `tautomer nH`: Tautomer NH ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_HBD_TAUTOMERIC_NH))
    - `oxygen amide term`: Oxygen amide terminal \
        ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_HBD_OXYGEN_AMIDE_TERM))
"""

AROMATIC_SMARTS = {
    "arom_4": SMARTS_STR_AROMATIC_4,
    "arom_5": SMARTS_STR_AROMATIC_5,
    "arom_6": SMARTS_STR_AROMATIC_6,
    "arom_7": SMARTS_STR_AROMATIC_7,
    "arom_8": SMARTS_STR_AROMATIC_8,
}
"""
A dictionary of aromatic SMARTS patterns. The number in the key indicates the number of atoms in the aromatic ring.
!!! tip "Definitions"
    - `arom_4`: Aromatic 4 ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_AROMATIC_4))
    - `arom_5`: Aromatic 5 ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_AROMATIC_5))
    - `arom_6`: Aromatic 6 ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_AROMATIC_6))
    - `arom_7`: Aromatic 7 ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_AROMATIC_7))
    - `arom_8`: Aromatic 8 ([source](smarts_defs.md#lahuta.config.smarts_patterns.SMARTS_STR_AROMATIC_8))
"""
