import sys
from _lahuta import Luni

luni = Luni(sys.argv[1])
ns = luni.get_neighbors()
print(ns.size())
print(len(ns.get_neighbors()))
print(ns.get_neighbors()[:10])
print(luni.get_cutoff())
