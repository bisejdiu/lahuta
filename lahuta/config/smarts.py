AMIDE_SMARTS = "[NX3][CX3](=[OX1])[#6]"  # DEFINITION FROM `http://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html`


ATOM_TYPES = {
    "hbond_acceptor": {
        "acceptor": "[#8,#9,$([#16;H0,H1;v2,v1]),$([N;v3;!$(N-*=!@[O,N,P,S]);!$(N-!@a);!$([NH]=!@*)]),$([nH0;+0])]",
        "enol": "[$([nH]:@c(=O))]",
        "tautomeric nH": "[$([n;H1;v3;!$([nH]cccc)])]",
        # AMBIGUITY OF TERMINAL AMIDES MAY AFFECT NON-PROTEIN AMIDES
        "NH2 terminal amide": "[$([N;H2;v3;$(N-C(=O))])]",
    },
    "hbond_donor": {
        "donor": "[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]",
        "oxygen acid": "[$([O;H0;$(O=C([OH])-*)])]",
        "tautomer nH": "[$(n:a:[nH])]",
        # AMBIGUITY OF TERMINAL AMIDES MAY AFFECT NON-PROTEIN AMIDES
        "oxygen amide term": "[$([O;H0;$(O=C-[NH2])])]",
    },
    "xbond_acceptor": {
        # SAME AS HBA
        "acceptor": "[#8,#9,$([#16;H0,H1;v2,v1]),$([N;v3;!$(N-*=!@[O,N,P,S]);!$(N-!@a);!$([NH]=!@*)]),$([nH0;+0])]",
        "enol": "[$([nH]:@c(=O))]",
        "tautomeric nH": "[$([n;H1;v3;!$([nH]cccc)])]",
        # AMBIGUITY OF TERMINAL AMIDES MAY AFFECT NON-PROTEIN AMIDES
        "NH2 terminal amide": "[$([N;H2;v3;$(N-C(=O))])]",
    },
    "xbond_donor": {"donor": "[Cl,Br,I;X1;$([Cl,Br,I]-[#6])]"},
    "weak_hbond_acceptor": {
        # SAME AS HBA
        "acceptor": "[#8,#9,$([#16;H0,H1;v2,v1]),$([N;v3;!$(N-*=!@[O,N,P,S]);!$(N-!@a);!$([NH]=!@*)]),$([nH0;+0])]",
        "enol": "[$([nH]:@c(=O))]",
        "tautomeric nH": "[$([n;H1;v3;!$([nH]cccc)])]",
        # AMBIGUITY OF TERMINAL AMIDES MAY AFFECT NON-PROTEIN AMIDES
        "NH2 terminal amide": "[$([N;H2;v3;$(N-C(=O))])]",
        "c-x halogens": "[Cl,Br,I;X1;$([Cl,Br,I]-[#6])]",
    },
    "weak_hbond_donor": {"donor": "[#6!H0]"},
    # SEE RDKIT `BaseFeatures.fdef`
    "pos_ionisable": {
        "rdkit basic group": "[$([N;H2&+0][C;!$(C=*)]),$([N;H1&+0]([C;!$(C=*)])[C;!$(C=*)]),$([N;H0&+0]([C;!$(C=*)])([C;!$(C=*)])[C;!$(C=*)]);!$(N[a])]",
        "imidazole": "[n;R1]1[c;R1][n;R1][c;R1][c;R1]1",
        "guanidine amidine": "NC(=N)",
        "rdkit posn": "[#7;+;!$([N+]-[O-])]",
        "cations": "[$([*+1,*+2,*+3]);!$([N+]-[O-])]",
        "metals": "[Li,Be,Na,Mg,Al,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Rb,Sr,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Cs,Ba,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi,Po,Fr,Ra,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf]",
    },
    "neg_ionisable": {
        "O acidic group": "[$([OH,O-]-[C,S,N,P,Cl,Br,I]=O),$(O=[C,S,N,P,Cl,Br,I]-[OH,O-])]",
        "anions": "[*-1,*-2]",
    },
    "hydrophobe": {"hydrophobe": "[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,Cl+0,Br+0,I+0]"},
    "carbonyl_oxygen": {"oxygen": "[$([OH0]=[CX3,c]);!$([OH0]=[CX3,c]-[OH,O-])]"},
    "carbonyl_carbon": {"carbon": "[$([CX3,c]=[OH0]);!$([CX3,c](=[OH0])-[OH,O-])]"},
    "aromatic": {
        "arom_4": "[a;r4,!R1&r3]1:[a;r4,!R1&r3]:[a;r4,!R1&r3]:[a;r4,!R1&r3]:1",
        "arom_5": "[a;r5,!R1&r4,!R1&r3]1:[a;r5,!R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:[a;r5,!R1&r4,!R1&r3]:1",
        "arom_6": "[a;r6,!R1&r5,!R1&r4,!R1&r3]1:[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:[a;r6,!R1&r5,!R1&r4,!R1&r3]:1",
        "arom_7": "[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]1:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:1",
        "arom_8": "[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]1:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:[a;r8,!R1&r7,!R1&r6,!R1&r5,!R1&r4,!R1&r3]:1",
    },
}
