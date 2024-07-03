import time
import numpy as np
from lahuta import Luni
from lahuta.utils.mda import mda_psuedobox_from_atomgroup
from MDAnalysis.lib.nsgrid import FastNS


def degsin(deg: float):
    deg = deg * np.pi / 180.0
    return np.sin(deg)


# file_name = "1kx2_small.cif"
# file_name = "8rt5_17k.cif"
# file_name = "8rdu_60k.cif"
# file_name = "8ugh_100k.cif"
# file_name = "8sf7_1M.cif"
file_name = "8glv_4M.cif"
luni = Luni("data/" + file_name)
print("NUMBER OF ATOMS: ", luni.n_atoms)
mda = luni.to("mda")

start = time.time()
positions, dimensions = mda_psuedobox_from_atomgroup(mda)
gs = FastNS(cutoff=5, coords=positions, box=dimensions, pbc=False)
r = gs.self_search()
end = time.time()
print(r.get_pairs().shape)
print(gs.get_counter())
print(end - start)
