# from .core.groups import AtomGroup

from lahuta import AtomGroup
from lahuta.core.universe import Universe
from lahuta.contacts.covalent import CovalentContactStrategy
from lahuta.contacts.metal import MetalContactStrategy
from lahuta.contacts.hbonds import HBondContactStrategy

if __name__ == "__main__":

    import numpy as np

    np.set_printoptions(suppress=True)

    u = Universe("/home/bisejdiu/p/lahuta/arpeggio/notebooks/1KX2.pdb")
    n = u.compute_neighbors()
    # print(n)
    m = (
        n.type_filter(["carbonyl oxygen"], 0)
        .type_filter(["carbonyl oxygen"], 1)
        .distance_filter(4.0)
    )
    # print(m)

    cc = CovalentContactStrategy(u, n)
    # print(cc.contacts("dataframe", "compact"))

    mc = MetalContactStrategy(u, n)
    # print(mc.contacts("dataframe", "compact"))

    # v = (
    #     n.type_filter("hbond donor", 0)
    #     .type_filter("hbond acceptor", 1)
    #     .contact_type("hbond")
    #     .hbond_distance_filter(col=1)
    #     .hbond_angle_filter(col=0)
    # )
    # print(v.result_array)

    hb = HBondContactStrategy(u, n)
    print(hb.pairs)
    print(hb.contacts("dataframe", "compact"))
