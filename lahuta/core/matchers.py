"""
This module contains classes for SMARTS pattern matching on molecules.

The SMARTS pattern matching classes are used to match SMARTS patterns to atoms in a molecule.
This is how we assign atom types to molecules.

Classes:
    SmartsMatcherBase: Abstract base class for SMARTS pattern matching.
    SmartsMatcher: Sequential SMARTS pattern matching.
    ParallelSmartsMatcher: Parallel SMARTS pattern matching.

"""
import os
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Dict, List, Tuple

import numpy as np
from numpy.typing import NDArray
from openbabel import openbabel as ob
from scipy.sparse import dok_matrix

from lahuta.config.atoms import STANDARD_AMINO_ACIDS
from lahuta.config.smarts import AVAILABLE_ATOM_TYPES as ATypes
from lahuta.config.smarts import SmartsPatternRegistry
from lahuta.lahuta_types.openbabel import MolType, ObSmartPatternType, OBSmartsPatternWrapper


class SmartsMatcherBase(ABC):
    """
    A base class for different implementations of SMARTS pattern matching on molecules.

    This abstract class needs to be inherited by any class that implements SMARTS pattern matching.
    The subclass must implement the compute method.
    """

    @abstractmethod
    def compute(self, mol: MolType) -> NDArray[np.int8]:
        """
        Abstract method for SMARTS pattern matching.

        Args:
            mol (MolType): A molecule object to match patterns on.

        Raises:
            NotImplementedError: This is an abstract method that needs to be implemented in the subclass.

        Returns:
            NDArray[np.int8]: An array of atom types as defined in the SmartsPatternRegistry dictionary.
        """
        raise NotImplementedError("Subclasses must implement this method")


class SmartsMatcher(SmartsMatcherBase):
    """
    Matches SMARTS patterns to atoms in a molecule.

    This class performs sequential SMARTS pattern matching on atoms in a molecule.
    It inherits from the SmartsMatcherBase abstract base class.
    """

    def compute(self, mol: MolType, mda) -> NDArray[np.int8]:
        """
        Performs SMARTS pattern matching on a molecule.

        Args:
            mol (MolType): A molecule object to match patterns on.

        Returns:
            NDArray[np.int8]: An array of atom types that match the SMARTS patterns in the given molecule.
        """

        shape = (mol.NumAtoms(), len(ATypes))
        shape = (mda.universe.atoms.n_atoms, len(ATypes))
        # atypes_array: NDArray[np.int8] = np.zeros(shape, dtype=np.int8)
        dok_atyps = dok_matrix(shape, dtype=np.int8)

        for atom_type in SmartsPatternRegistry:
            smartsdict = SmartsPatternRegistry[atom_type.name].value
            for smarts in smartsdict.values():
                ob_smart: ObSmartPatternType = OBSmartsPatternWrapper(ob.OBSmartsPattern())
                ob_smart.Init(str(smarts))
                ob_smart.Match(mol)

                matches = [x[0] for x in ob_smart.GetMapList()]
                for match in matches:
                    atom = mol.GetAtom(match)

                    if atom.GetResidue().GetName() not in STANDARD_AMINO_ACIDS:
                        if atom.GetId() in [1174, 1175, 1179, 1181]:
                            print('atom.GetId()', atom.GetId(), atom_type.name)
                        # atypes_array[atom.GetId(), ATypes[atom_type.name]] = 1
                        dok_atyps[atom.GetId(), ATypes[atom_type.name]] = 1

        return dok_atyps


class ParallelSmartsMatcher(SmartsMatcherBase):
    """
    Matches SMARTS patterns to atoms in a molecule using multiple threads.

    This class performs SMARTS pattern matching on atoms in a molecule using multiple threads for
    improved performance. It inherits from the SmartsMatcherBase abstract base class.
    """

    def __init__(self) -> None:
        self.precomputed_ob_smarts = self.precompute_ob_smarts()

    def precompute_ob_smarts(self) -> Dict[str, List[ObSmartPatternType]]:
        """
        Precomputes and stores the Open Babel SMARTS patterns for all atom types.

        Returns:
            Dict[str, List[ObSmartPatternType]]: A dictionary with atom type names as keys and lists of
                                                precomputed Open Babel SMARTS patterns as values.
        """

        precomputed_ob_smarts: Dict[str, List[ObSmartPatternType]] = {}
        for atom_type in SmartsPatternRegistry:
            smartsdict = SmartsPatternRegistry[atom_type.name].value
            precomputed_ob_smarts[atom_type.name] = []
            for smarts in smartsdict.values():
                ob_smart: ObSmartPatternType = OBSmartsPatternWrapper(ob.OBSmartsPattern())
                ob_smart.Init(str(smarts))
                precomputed_ob_smarts[atom_type.name].append(ob_smart)
        return precomputed_ob_smarts

    def match_ob_smarts(
        self,
        ob_smart: ObSmartPatternType,
        mol: MolType,
        atypes: Dict[str, int],
        atom_type: str,
    ) -> List[Tuple[Any, int]]:
        """
        Matches an Open Babel SMARTS pattern to a molecule.

        Args:
            ob_smart (ObSmartPatternType): An Open Babel SMARTS pattern.
            mol (MolType): A molecule object to match the pattern on.
            atypes (Dict[str, int]): A dictionary of atom types.
            atom_type (str): The name of the atom type that the SMARTS pattern represents.

        Returns:
            List[Tuple[Any, int]]: A list of tuples, where each tuple contains the matched atom's
                                    index and the corresponding atom type.
        """

        ob_smart.Match(mol)
        matches = [x[0] for x in ob_smart.GetMapList()]
        return [(match, atypes[atom_type]) for match in matches]

    def compute(self, mol: MolType) -> NDArray[np.int8]:
        """
        Performs SMARTS pattern matching on a molecule using multiple threads.

        Args:
            mol (MolType): A molecule object to match patterns on.

        Returns:
            NDArray[np.int8]: An array of atom types that match the SMARTS patterns in the given molecule.
        """

        shape = (mol.NumAtoms(), len(ATypes))
        atypes_array: NDArray[np.int8] = np.zeros(shape, dtype=np.int8)

        num_threads = os.cpu_count()

        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            for atom_type, ob_smarts_list in self.precomputed_ob_smarts.items():
                future_matches = [
                    executor.submit(self.match_ob_smarts, ob_smart, mol, ATypes, atom_type)
                    for ob_smart in ob_smarts_list
                ]

                for future in future_matches:
                    matches = future.result()
                    for match, atype in matches:
                        atom = mol.GetAtom(match)
                        if atom.GetResidue().GetName() not in STANDARD_AMINO_ACIDS:
                            atypes_array[atom.GetIdx() - 1, atype] = 1

        return atypes_array
