import os
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from openbabel import openbabel as ob

from lahuta.config.atoms import STANDARD_AMINO_ACIDS
from lahuta.config.smarts import AVAILABLE_ATOM_TYPES as ATypes


class SmartsMatcherBase(ABC):
    """
    Base class for SMARTS pattern matching on molecules.

    This class serves as an abstract base class for concrete implementations
    of SMARTS pattern matching, such as SmartsMatcher and ParallelSmartsMatcher.
    """

    def __init__(self, atom_types):
        self.atom_types = atom_types

    @abstractmethod
    def compute(self, mol):
        raise NotImplementedError("Subclasses must implement this method")


class SmartsMatcher(SmartsMatcherBase):
    """
    A class for matching SMARTS patterns to atoms in a molecule.

    This class performs sequential SMARTS pattern matching on atoms in a
    molecule, returning atom types as defined in the atom_types dictionary.
    Inherits from the SmartsMatcherBase abstract base class.
    """

    def compute(self, mol):
        atypes_array = np.zeros((mol.NumAtoms(), len(ATypes)))

        for atom_type in self.atom_types:
            smartsdict = self.atom_types[atom_type.name].value
            for smarts in smartsdict.values():
                ob_smart = ob.OBSmartsPattern()
                ob_smart.Init(str(smarts))
                ob_smart.Match(mol)

                matches = [x[0] for x in ob_smart.GetMapList()]
                for match in matches:
                    atom = mol.GetAtom(match)

                    if atom.GetResidue().GetName() not in STANDARD_AMINO_ACIDS:
                        atypes_array[atom.GetId(), ATypes[atom_type.name].value] = 1

        return atypes_array


class ParallelSmartsMatcher(SmartsMatcherBase):
    """
    A class for parallel matching of SMARTS patterns to atoms in a molecule.

    This class utilizes multiple threads for parallel SMARTS pattern matching
    on atoms in a molecule, returning atom types as defined in the atom_types
    dictionary. Inherits from the SmartsMatcherBase abstract base class.
    """

    def __init__(self, atom_types, n_atoms):
        super().__init__(atom_types, n_atoms)
        self.precomputed_ob_smarts = self.precompute_ob_smarts()

    def precompute_ob_smarts(self):
        precomputed_ob_smarts = {}
        for atom_type in self.atom_types:
            smartsdict = self.atom_types[atom_type.name].value
            precomputed_ob_smarts[atom_type] = []
            for smarts in smartsdict.values():
                ob_smart = ob.OBSmartsPattern()
                ob_smart.Init(str(smarts))
                precomputed_ob_smarts[atom_type].append(ob_smart)
        return precomputed_ob_smarts

    def match_ob_smarts(self, ob_smart, mol, atypes, atom_type):
        ob_smart.Match(mol)
        matches = [x[0] for x in ob_smart.GetMapList()]
        return [(match, atypes[atom_type.name].value) for match in matches]

    def compute(self, mol):
        atypes_array = np.zeros((mol.NumAtoms(), len(ATypes)))

        num_threads = os.cpu_count()

        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            for atom_type, ob_smarts_list in self.precomputed_ob_smarts.items():
                future_matches = [
                    executor.submit(
                        self.match_ob_smarts, ob_smart, mol, ATypes, atom_type
                    )
                    for ob_smart in ob_smarts_list
                ]

                for future in future_matches:
                    matches = future.result()
                    for match, atype in matches:
                        atom = mol.GetAtom(match)
                        if atom.GetResidue().GetName() not in STANDARD_AMINO_ACIDS:
                            atypes_array[atom.GetIdx() - 1, atype] = 1

        return atypes_array
