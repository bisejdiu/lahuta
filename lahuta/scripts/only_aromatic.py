import time

import MDAnalysis as mda

from lahuta.contacts import AtomPlaneContacts, F
from lahuta.contacts.plane_plane import PlanePlaneContacts  # APDataFrameFactory,; PPDataFrameFactory,
from lahuta.core.universe import Universe

if __name__ == "__main__":
    # Load the universe

    # u = Universe("/home/bisejdiu/2023/lahuta/lahuta/tests/data/1KX2.pdb")
    # u = Universe("/home/bisejdiu/tutorials/lahuta-notebooks/data/1KX2_rcsb.pdb")
    u = Universe("/home/bisejdiu/tutorials/lahuta-notebooks/data/1KX2_rcsb.cif")
    # mda_u = mda.Universe("/home/bisejdiu/tutorials/lahuta-notebooks/data/1KX2_rcsb.pdb")
    # # mda_u_pp = mda_u.select_atoms("(protein and not resname ARG) or resname HEC")
    # # mda_u_pp = mda_u.select_atoms("protein and not resname ARG")
    # mda_u_pp = mda_u.select_atoms("all")
    # u = Universe(mda_u.atoms)
    # u = Universe("/home/bisejdiu/tutorials/lahuta-notebooks/data/4GSW.pdb")
    # u = Universe("/home/bisejdiu/tutorials/lahuta-notebooks/data/4GSW.cif")
    # end = time.time()
    # u = Universe("/home/bisejdiu/tutorials/lahuta-notebooks/data/8djb.cif")
    start = time.time()
    n = u.compute_neighbors(res_dif=2)
    print("Finished computing neighbors", n.pairs.shape)

    # ionic = F.ionic_neighbors(n)
    # print(ionic.pairs.shape, "ionic")

    aromatic = F.aromatic_neighbors(n)
    print(aromatic.pairs.shape, "aromatic")

    end = time.time()
    print("Time elapsed: ", end - start)
