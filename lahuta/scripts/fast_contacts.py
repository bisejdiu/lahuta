import time

import lahuta.contacts.contacts as contacts
from lahuta.contacts.plane import (APDataFrameFactory, AtomPlaneContacts,
                                   PlanePlaneContacts, PPDataFrameFactory)
from lahuta.core.universe import Universe

if __name__ == "__main__":
    # Load the universe
    u = Universe("/home/bisejdiu/2023/lahuta/lahuta/tests/data/1KX2.pdb")
    # u = Universe("/home/bisejdiu/tutorials/lahuta-notebooks/data/1cbs.cif")
    n = u.compute_neighbors()

    start = time.time()
    # Compute contacts
    cov = contacts.covalent_neighbors(n)
    print(cov.pairs.shape)

    met = contacts.metalic_neighbors(n)
    print(met.pairs.shape)

    carb = contacts.carbonyl_neighbors(n)
    print(carb.pairs.shape)

    hb = contacts.hbond_neighbors(n)
    print(hb.pairs.shape)

    whb = contacts.weak_hbond_neighbors(n)
    print(whb.pairs.shape)

    ionic = contacts.ionic_neighbors(n)
    print(ionic.pairs.shape)

    aromatic = contacts.aromatic_neighbors(n)
    print(aromatic.pairs.shape)

    hydrophobic = contacts.hydrophobic_neighbors(n)
    print(hydrophobic.pairs.shape)

    phb = contacts.polar_hbond_neighbors(n)
    print(phb.pairs.shape)

    wphb = contacts.weak_polar_hbond_neighbors(n)
    print(wphb.pairs.shape)

    vdw = contacts.vdw_neighbors(n)
    print(vdw.pairs.shape)

    ap = AtomPlaneContacts(u)
    ap.compute_contacts()

    # z = ap.get_contacts()
    # print(z)
    cp = ap.carbon_pi.contacts(ap.neighbors, ap.angles)
    print(cp.pairs.shape)

    cp2 = ap.cation_pi.contacts(ap.neighbors, ap.angles)
    print(cp2.pairs.shape)

    dp = ap.donor_pi.contacts(ap.neighbors, ap.angles)
    print(dp.pairs.shape)

    sp = ap.sulphur_pi.contacts(ap.neighbors, ap.angles)
    print(sp.pairs.shape)

    pp = PlanePlaneContacts(u)
    pp.compute_contacts()
    print(pp.pairs.shape)

    print(PPDataFrameFactory(pp, df_format="expanded").dataframe())

    end = time.time()
    print("Time elapsed: ", end - start)
