# from .core.groups import AtomGroup

from lahuta import AtomGroup
from lahuta.core.universe import Universe
from lahuta.contacts.covalent import CovalentContactStrategy
from lahuta.contacts.metal import MetalContactStrategy
from lahuta.contacts.hbonds import HBondContactStrategy
from lahuta.contacts.whbond import WeakHBondContactStrategy
from lahuta.contacts.ionic import IonicContactStrategy
from lahuta.config import config


if __name__ == "__main__":

    import numpy as np

    np.set_printoptions(suppress=True)

    u = Universe("/home/bisejdiu/p/lahuta/lahuta/notebooks/1KX2.pdb")
    n = u.compute_neighbors()
    # print(n)
    # m = (
    #     n.type_filter(["carbonyl oxygen"], 0)
    #     .type_filter(["carbonyl oxygen"], 1)
    #     .distance_filter(4.0)
    # )
    # print(m)

    # cc = CovalentContactStrategy(u, n)
    # print(cc.contacts("dataframe", "compact"))

    # mc = MetalContactStrategy(u, n)
    # print(mc.contacts("dataframe", "compact"))

    # v = (
    #     n.type_filter("hbond donor", 0)
    #     .type_filter("hbond acceptor", 1)
    #     .contact_type("hbond")
    #     .hbond_distance_filter(col=1)  # type: ignore
    #     .hbond_angle_filter(col=0)  # type: ignore
    # )
    # print(v.result_array)

    # hb = HBondContactStrategy(u, n)
    # whb = WeakHBondContactStrategy(u, n)
    # # print(hb.pairs)
    # print(whb.contacts("dataframe", "print").head(40))
    # print(whb.contacts("dataframe", "print").shape)

    # v1 = (
    #     n.type_filter("pos ionisable", 0)
    #     .type_filter("neg ionisable", 1)
    #     .distance_filter(config.CONTACT_TYPES["ionic"]["distance"])
    # )
    # v2 = (
    #     n.type_filter("neg ionisable", 0)
    #     .type_filter("pos ionisable", 1)
    #     .distance_filter(config.CONTACT_TYPES["ionic"]["distance"])
    # )

    # print(v1.pairs)
    # print(v2.pairs)

    i = IonicContactStrategy(u, n)
    print(i.contacts("dataframe", "print").head(40))
    print(i.contacts("dataframe", "print").shape)
