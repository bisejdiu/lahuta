from abc import ABC, abstractmethod

import numpy as np

from lahuta.config.atoms import PROT_ATOM_TYPES


class ProteinTypeAssignerBase(ABC):
    """
    Base class for assigning protein atom types.

    This class serves as an abstract base class for concrete implementations
    of protein atom type assignment, such as VectorizedProteinTypeAssigner
    and LegacyProteinTypeAssigner.
    """

    def __init__(self, protein_atomgroup, ta):
        self.protein_atomgroup = protein_atomgroup
        self.ta = ta

    @abstractmethod
    def compute(self, atypes_array, atypes):
        """
        Abstract method for computing protein atom types.
        """
        raise NotImplementedError


class VectorizedProteinTypeAssigner(ProteinTypeAssignerBase):
    """
    A class for vectorized assignment of protein atom types.

    This class utilizes efficient NumPy array manipulation to assign atom types
    to proteins in a vectorized manner. Inherits from the ProteinTypeAssignerBase
    abstract base class.
    """

    def compute(self, atypes_array, atypes):
        # print("atypes", atypes)
        # resname, atom_name = self.ta["resname"], self.ta["name"]

        # resname_str = resname[self.protein_atomgroup.resindices].astype(str)
        # atom_name_str = atom_name[self.protein_atomgroup.indices].astype(str)

        resname_str = self.protein_atomgroup.resnames.astype(str)
        atom_name_str = self.protein_atomgroup.names.astype(str)

        atom_id_labels = np.core.defchararray.add(  # type: ignore
            np.core.defchararray.strip(resname_str),  # type: ignore
            np.core.defchararray.strip(atom_name_str),  # type: ignore
        )

        # prot_atom_types_array = [
        #     list(atom_ids_label_set) for atom_ids_label_set in PROT_ATOM_TYPES.values()
        # ]
        prot_atom_types_array = [list(PROT_ATOM_TYPES[key]) for key in atypes.keys()]
        # prot_atom_types_array = []
        # for key in atypes.keys():
        #     print("adding", key, atypes[key])
        #     prot_atom_types_array.append(list(PROT_ATOM_TYPES[key]))

        mask = np.array(
            [
                np.isin(atom_id_labels, prot_atom_types)
                for prot_atom_types in prot_atom_types_array
            ]
        )

        true_indices = np.argwhere(mask)

        atypes_array[true_indices[:, 1], true_indices[:, 0]] = 1

        return atypes_array


class LegacyProteinTypeAssigner(ProteinTypeAssignerBase):
    """
    A class for legacy assignment of protein atom types.

    This class uses a less efficient, loop-based approach to assign atom types
    to proteins. Inherits from the ProteinTypeAssignerBase abstract base class.
    """

    def compute(self, atypes_array, atypes):
        for residue in self.protein_atomgroup.residues:
            for atom in residue.atoms:
                for atom_type in list(PROT_ATOM_TYPES.keys()):
                    atypes_array[atom.index, atypes[atom_type]] = 0

                for atom_type, atom_ids in PROT_ATOM_TYPES.items():
                    atom_id = residue.resname.strip() + atom.name.strip()
                    if atom_id in atom_ids:
                        atypes_array[atom.index, atypes[atom_type]] = 1

        return atypes_array
