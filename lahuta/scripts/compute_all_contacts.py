import time

from lahuta import AtomGroup
from lahuta.config.defaults import CONTACTS
from lahuta.contacts import (AromaticContacts, CarbonylContacts,
                             CovalentContacts, HBondContacts,
                             HydrophobicContacts, IonicContacts, MetalContacts,
                             PolarHBondContacts, WeakHBondContacts,
                             WeakPolarHBondContacts)
from lahuta.contacts.plane import (APDataFrameFactory, AtomPlaneContacts,
                                   EnumContactStrategies, PlanePlaneContacts,
                                   PPDataFrameFactory)
from lahuta.contacts.vdw import VanDerWaalsContacts
from lahuta.core.universe import Universe

# measure time
start = time.time()
# Load the universe
u = Universe("/home/bisejdiu/p/lahuta/lahuta/notebooks/2RH1.pdb")
n = u.compute_neighbors()

# Compute contacts
cc = CovalentContacts(u, n)
# c = cc.contacts("dataframe", "expanded")
# print("covalent", c.shape)
# print(c)

mc = MetalContacts(u, n)
# c = mc.contacts("dataframe", "expanded")
# print("metal", c.shape)

cc = CarbonylContacts(u, n)
# c = cc.contacts("dataframe", "expanded")
# print("carbonyl", c.shape)

hb = HBondContacts(u, n)
# c = hb.contacts("dataframe", "expanded")
# print("hbond", c.shape)

whb = WeakHBondContacts(u, n)
# c = whb.contacts("dataframe", "expanded")
# print("weak hbond", c.shape)

i = IonicContacts(u, n)
# c = i.contacts("dataframe", "expanded")
# print("ionic", c.shape)

a = AromaticContacts(u, n)
# c = a.contacts("dataframe", "expanded")
# print("aromatic", c.shape)

h = HydrophobicContacts(u, n)
# c = h.contacts("dataframe", "expanded")
# print("hydrophobic", c.shape)

ph = PolarHBondContacts(u, n)
# c = ph.contacts("dataframe", "expanded")
# print("polar hbond", c.shape)

wph = WeakPolarHBondContacts(u, n)
# c = ph.contacts("dataframe", "expanded")
# print("weak polar hbond", c.shape)

vdw = VanDerWaalsContacts(u, n)
# c = vdw.contacts("dataframe", "expanded")
# print("vdw", c.shape)

end = time.time()
print("time", end - start)


ap = AtomPlaneContacts(u)
ap.compute_contacts()

s = ap.call_strategy(EnumContactStrategies.DonorPI_contacts)
c = APDataFrameFactory(
    s._ua.atoms[s.pairs], s.distances, ap.rings, "DonorPI", "expanded"
).dataframe()
print("donor pi", c.shape)

s = ap.call_strategy(EnumContactStrategies.CarbonPI_contacts)
c = APDataFrameFactory(
    s._ua.atoms[s.pairs], s.distances, ap.rings, "CarbonPI", "expanded"
).dataframe()
print("carbon pi", c.shape)

s = ap.call_strategy(EnumContactStrategies.CationPI_contacts)
c = APDataFrameFactory(
    s._ua.atoms[s.pairs], s.distances, ap.rings, "CationPI", "expanded"
).dataframe()
print("cation pi", c.shape)

s = ap.call_strategy(EnumContactStrategies.SulphurPI_contacts)
c = APDataFrameFactory(
    s._ua.atoms[s.pairs], s.distances, ap.rings, "SulphurPI", "expanded"
).dataframe()
print("sulphur pi", c.shape)

pp = PlanePlaneContacts(u)
_ = pp.compute_contacts()

c = PPDataFrameFactory(
    pp.ua.atoms[pp.pairs], pp.distances, pp.contact_labels, "expanded"
).dataframe()
print("plane plane", c.shape)

# end = time.time()
# print("time", end - start)
