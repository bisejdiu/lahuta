import time

import MDAnalysis as mda

import lahuta.contacts.contacts as contacts
from lahuta.contacts.plane import (
    APDataFrameFactory,
    AtomPlaneContacts,
    PlanePlaneContacts,
    PPDataFrameFactory,
)
from lahuta.core.universe import Universe

if __name__ == "__main__":
    # Load the universe

    # u = Universe("/home/bisejdiu/2023/lahuta/lahuta/tests/data/1KX2.pdb")
    # u = Universe("/home/bisejdiu/tutorials/lahuta-notebooks/data/1KX2_rcsb.pdb")
    # u = Universe("/home/bisejdiu/tutorials/lahuta-notebooks/data/1KX2_rcsb.cif")
    mda_u = mda.Universe("/home/bisejdiu/tutorials/lahuta-notebooks/data/1KX2_rcsb.pdb")
    # mda_u_pp = mda_u.select_atoms("(protein and not resname ARG) or resname HEC")
    mda_u_pp = mda_u.select_atoms("protein and not resname ARG")
    # mda_u_pp = mda_u.select_atoms("all")
    u = Universe(mda_u_pp.atoms)
    # u = Universe("/home/bisejdiu/tutorials/lahuta-notebooks/data/4GSW.pdb")
    # u = Universe("/home/bisejdiu/tutorials/lahuta-notebooks/data/4GSW.cif")
    # end = time.time()
    # u = Universe("/home/bisejdiu/tutorials/lahuta-notebooks/data/8djb.cif")
    start = time.time()
    n = u.compute_neighbors(res_dif=2)
    print("Finished computing neighbors")

    # Compute contacAMIDE_SMARTSts
    cov = contacts.covalent_neighbors(n)
    print(cov.pairs.shape, "cov")

    met = contacts.metalic_neighbors(n)
    print(met.pairs.shape, "met")

    carb = contacts.carbonyl_neighbors(n)
    print(carb.pairs.shape, "carb")

    hb = contacts.hbond_neighbors(n)
    print(hb.pairs.shape, "hb")

    whb = contacts.weak_hbond_neighbors(n)
    print(whb.pairs.shape, "whb")

    ionic = contacts.ionic_neighbors(n)
    print(ionic.pairs.shape, "ionic")

    aromatic = contacts.aromatic_neighbors(n)
    print(aromatic.pairs.shape, "aromatic")

    hydrophobic = contacts.hydrophobic_neighbors(n)
    print(hydrophobic.pairs.shape, "hydrophobic")

    phb = contacts.polar_hbond_neighbors(n)
    print(phb.pairs.shape, "phb")

    wphb = contacts.weak_polar_hbond_neighbors(n)
    print(wphb.pairs.shape, "wphb")

    vdw = contacts.vdw_neighbors(n)
    print(vdw.pairs.shape, "vdw")

    ap = AtomPlaneContacts(u)
    ap.compute_contacts()

    cp = ap.carbon_pi.contacts(ap.neighbors, ap.angles)  # type: ignore
    print(cp.pairs.shape, "cp")

    cp2 = ap.cation_pi.contacts(ap.neighbors, ap.angles)  # type: ignore
    print(cp2.pairs.shape, "cp2")

    dp = ap.donor_pi.contacts(ap.neighbors, ap.angles)  # type: ignore
    print(dp.pairs.shape, "dp")

    sp = ap.sulphur_pi.contacts(ap.neighbors, ap.angles)  # type: ignore
    print(sp.pairs.shape, "sp")

    pp = PlanePlaneContacts(u)
    pp.compute_contacts()
    print(pp.pairs.shape, "pp")

    print(PPDataFrameFactory(pp, df_format="expanded").dataframe())

    end = time.time()
    print("Time elapsed: ", end - start)
