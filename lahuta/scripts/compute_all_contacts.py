import time

from plane_plane import APDataFrameFactory, AtomPlaneContacts, PlanePlaneContacts, PPDataFrameFactory

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
from lahuta.contacts.vdw import VanDerWaalsContacts
from lahuta.core.universe import Luni

# measure time
start = time.time()
# Load the universe
u = Luni("/home/bisejdiu/p/lahuta/lahuta/notebooks/2RH1.pdb")
n = u.compute_neighbors()

# Compute contacts
cc = CovalentContacts(u, n)
c = cc.contacts("dataframe", "expanded")
print("covalent", c.shape)
print(c)

mc = MetalContacts(u, n)
c = mc.contacts("dataframe", "expanded")
print("metal", c.shape)

cc = CarbonylContacts(u, n)
c = cc.contacts("dataframe", "expanded")
print("carbonyl", c.shape)

hb = HBondContacts(u, n)
c = hb.contacts("dataframe", "expanded")
print("hbond", c.shape)

whb = WeakHBondContacts(u, n)
c = whb.contacts("dataframe", "expanded")
print("weak hbond", c.shape)

i = IonicContacts(u, n)
c = i.contacts("dataframe", "expanded")
print("ionic", c.shape)

a = AromaticContacts(u, n)
c = a.contacts("dataframe", "expanded")
print("aromatic", c.shape)

h = HydrophobicContacts(u, n)
c = h.contacts("dataframe", "expanded")
print("hydrophobic", c.shape)

ph = PolarHBondContacts(u, n)
c = ph.contacts("dataframe", "expanded")
print("polar hbond", c.shape)

wph = WeakPolarHBondContacts(u, n)
c = ph.contacts("dataframe", "expanded")
print("weak polar hbond", c.shape)

vdw = VanDerWaalsContacts(u, n)
c = vdw.contacts("dataframe", "expanded")
print("vdw", c.shape)


end = time.time()
print("time", end - start)
