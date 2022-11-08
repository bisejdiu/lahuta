# from .core.groups import AtomGroup

from lahuta import AtomGroup
from lahuta.core.universe import Universe
from lahuta.contacts.covalent import CovalentContacts
from lahuta.contacts.metal import MetalContacts
from lahuta.contacts.hbonds import HBondContacts
from lahuta.contacts.whbond import WeakHBondContacts
from lahuta.contacts.ionic import IonicContacts
from lahuta.contacts.carbonyl import CarbonylContacts
from lahuta.contacts.aromatic import AromaticContacts
from lahuta.contacts.hydrophobic import HydrophobicContacts
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

    # cc = CovalentContacts(u, n)
    # print(cc.contacts("dataframe", "compact"))

    # mc = MetalContacts(u, n)
    # print(mc.contacts("dataframe", "compact"))

    # cc = CarbonylContacts(u, n)
    # # print(cc.pairs)
    # print(cc.contacts("dataframe", "compact"))

    # v = (
    #     n.type_filter("hbond donor", 0)
    #     .type_filter("hbond acceptor", 1)
    #     .contact_type("hbond")
    #     .hbond_distance_filter(col=1)  # type: ignore
    #     .hbond_angle_filter(col=0)  # type: ignore
    # )
    # print(v.result_array)

    # hb = HBondContacts(u, n)
    # # whb = WeakHBondContacts(u, n)
    # # print(hb.pairs)
    # # ww = hb.contacts("dataframe", "print")
    # # ww.head(20)
    # print(hb.contacts("dataframe", "print").head(40))
    # print(hb.contacts("dataframe", "print").shape)

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

    # i = IonicContacts(u, n)
    # print(i.col1.indices, "indices")
    # print(i.contacts("dataframe", "print").head(40))
    # print(i.contacts("dataframe", "print").shape)

    # a = AromaticContacts(u, n)
    # print(a.contacts("dataframe", "print").head(40))
    # print(a.contacts("dataframe", "print").shape)

    h = HydrophobicContacts(u, n)
    print(h.contacts("dataframe", "print").head(40))
    print(h.contacts("dataframe", "print").shape)
