# from .core.groups import AtomGroup

from lahuta import AtomGroup
from lahuta.config.defaults import CONTACTS
from lahuta.contacts import (
    AromaticContacts,
    CarbonylContacts,
    CovalentContacts,
    HBondContacts,
    HydrophobicContacts,
    IonicContacts,
    MetalContacts,
    PolarHBondContacts,
    WeakHBondContacts,
    WeakPolarHBondContacts,
)
from lahuta.contacts.atomplane import AtomPlaneContacts
from lahuta.contacts.planeplane import PlanePlaneContacts
from lahuta.contacts.vdw import VanDerWaalsContacts
from lahuta.core.universe import Universe

if __name__ == "__main__":
    import numpy as np

    np.set_printoptions(suppress=True)

    u = Universe("/home/bisejdiu/p/lahuta/lahuta/notebooks/1KX2.pdb")
    n = u.compute_neighbors()
    # print(n)
    # m = (
    #     n.type_filter(["carbonyl_oxygen"], 0)
    #     .type_filter(["carbonyl_oxygen"], 1)
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
    #     n.type_filter("hbond_donor", 0)
    #     .type_filter("hbond_acceptor", 1)
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
    #     n.type_filter("pos_ionisable", 0)
    #     .type_filter("neg_ionisable", 1)
    #     .distance_filter(CONTACTS["ionic"]["distance"])
    # )
    # v2 = (
    #     n.type_filter("neg_ionisable", 0)
    #     .type_filter("pos_ionisable", 1)
    #     .distance_filter(CONTACTS["ionic"]["distance"])
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

    # h = HydrophobicContacts(u, n)
    # print(h.contacts("dataframe", "print").head(40))
    # print(h.contacts("dataframe", "print").shape)

    # ph = PolarHBondContacts(u, n)
    # # print(ph.pairs)
    # pairs = ph.pairs
    # new_sorted_pairs = []
    # for pair in pairs:
    #     key = tuple(sorted(pair))
    #     new_sorted_pairs.append(key)

    # # sort according to the first element of the pair
    # new_sorted_pairs = sorted(new_sorted_pairs, key=lambda x: x[0])
    # print(np.array(new_sorted_pairs))
    # # print(ph.contacts("dataframe", "print").head(40))
    # print(ph.contacts("dataframe", "print").shape)

    # ph = WeakPolarHBondContacts(u, n)
    # pairs = ph.pairs
    # new_sorted_pairs = []
    # for pair in pairs:
    #     key = tuple(sorted(pair))
    #     new_sorted_pairs.append(key)

    # new_sorted_pairs = sorted(new_sorted_pairs, key=lambda x: x[0])
    # print(np.array(new_sorted_pairs))
    # print(ph.contacts("dataframe", "print"))
    # print(ph.contacts("dataframe", "compact"))
    # print(ph.contacts("dataframe", "expanded"))

    # vdw = VanDerWaalsContacts(u, n)
    # vdw.pairs = vdw.compute_contacts(remove_clashes=False)
    # vdw.distances = vdw.filter_distances()
    # # print(vdw.pairs)
    # print(vdw.contacts("dataframe", "compact"))

    # ap = AtomPlaneContacts(u)
    # z = ap.compute_contacts()
    # print(z)
    # print(len(z))

    # pp = PlanePlaneContacts(u)
    # z = pp.compute_contacts()
    # print(len(z))
    # print(z[0])
    # print(z[1])
