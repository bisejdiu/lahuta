import time

from lahuta.contacts import F
from lahuta.core.universe import Universe

if __name__ == "__main__":
    # Load the universe
    u = Universe("/home/bisejdiu/tutorials/lahuta-notebooks/data/1KX2_rcsb.cif")
    start = time.time()
    n = u.compute_neighbors(res_dif=2)
    print("Finished computing neighbors", n.pairs.shape)

    # ionic = F.ionic_neighbors(n)
    # print(ionic.pairs.shape, "ionic")

    aromatic = F.aromatic_neighbors(n)
    print(aromatic.pairs.shape, "aromatic")

    end = time.time()
    print("Time elapsed: ", end - start)
