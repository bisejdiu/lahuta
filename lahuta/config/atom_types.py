from collections import OrderedDict

AVAILABLE_ATOM_TYPES = OrderedDict(
    {
        "hbond_acceptor": 0,
        "pos_ionisable": 1,
        "carbonyl_oxygen": 2,
        "weak_hbond_donor": 3,
        "carbonyl_carbon": 4,
        "weak_hbond_acceptor": 5,
        "hbond_donor": 6,
        "neg_ionisable": 7,
        "aromatic": 8,
        "xbond_acceptor": 9,
        "hydrophobe": 10,
    }
)
