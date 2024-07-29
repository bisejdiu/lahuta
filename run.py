import sys
from _lahuta import Luni

luni = Luni(sys.argv[1])
ns = luni.getNeighborResults()
print(ns.getNeighborPairsSize())
print(len(ns.getNeighbors()))
print(ns.getNeighbors()[:10])
print(luni.getCutoff())
