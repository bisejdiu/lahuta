import MDAnalysis as mda
import numpy as np
import pandas as pd

from lahuta.core._loaders import GemmiLoader, TopologyLoader
from lahuta.core.base import FileLoader
from lahuta.core.obmol import OBMol


class StructureLoader:
    def __init__(self, file_path, is_pdb=False):
        self.file_path = file_path
        # extension = file_path.split(".")[-1]
        self.loader = GemmiLoader(file_path, is_pdb=is_pdb)
        # chains, residues, atoms = self.loader.create()

        assert self.loader.chains is not None
        assert self.loader.residues is not None
        assert self.loader.atoms is not None

    def load(self):
        mol = self.loader.to("mol")
        universe = self.loader.to("mda")
        return mol, universe


class PDBLoader(FileLoader):
    def load(self, *args):
        self._load_obabel()
        universe = mda.Universe(self.file_name, *args)
        return self.mol, universe
