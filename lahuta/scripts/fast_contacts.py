import time

import lahuta.contacts.contacts as contacts
from lahuta.core.universe import Universe

if __name__ == "__main__":
    # Load the universe
    u = Universe("/home/bisejdiu/p/lahuta/lahuta/notebooks/2RH1.pdb")
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

    end = time.time()
    print("Time elapsed: ", end - start)
