from typing import Dict

# http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
VDW_RADII = {"H": 1.2}

CONTACTS_DIST_MAX = 4.5

THETA_REQUIRED = set(["CARBONPI", "CATIONPI", "DONORPI", "HALOGENPI"])

CONTACTS: Dict[str, Dict[str, float]] = {
    "hbond": {
        "distance": 3.9,
        "polar distance": 3.5,
        "angle rad": 1.57,
        "angle degree": 90.0,
    },
    "weak hbond": {
        "distance": 3.6,
        "weak polar distance": 3.5,
        "angle rad": 2.27,
        "angle degree": 130.0,
        "cx angle min rad": 0.52,
        "cx angle min degree": 30.0,
        "cx angle max rad": 2.62,
        "cx angle max degree": 150.0,
    },
    "aromatic": {
        "distance": 4.0,
        "centroid_distance": 6.0,
        "atom_aromatic_distance": 4.5,
        "met_sulphur_aromatic_distance": 6.0,
    },
    "amide": {"centroid_distance": 6.0, "angle degree": 30.0, "angle rad": 0.52},
    "xbond": {
        "catmap distance": 1.85,  # SAME AS BROMINE
        "angle theta 1 rad": 2.09,
        "angle theta 1 degree": 120.0,
        "angle theta 2 min rad": 1.22,
        "angle theta 2 max rad": 2.97,
        "angle theta 2 min degree": 70.0,
        "angle theta 2 max degree": 170.0,
    },
    "ionic": {"distance": 4.0},
    "hydrophobic": {"distance": 4.5},
    "carbonyl": {"distance": 3.6},
    "metal": {"distance": 2.8},
}

GEMMI_SUPPRTED_FORMATS = {"cif", "mmcif", "cif.gz", "pdb", "pdb.gz"}
