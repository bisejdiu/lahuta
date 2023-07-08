import numpy as np


class Residues:
    def __init__(self, block_data):
        self.block_data = block_data

        self._resnames = np.array([], dtype="<U10")
        self._resids = np.array([], dtype=int)
        self._resindices = np.array([], dtype=int)
        self._icodes = np.array([], dtype="<U10")

        self._process()

    def _process(self):
        self._resnames = np.array(self.block_data.get("label_comp_id"))
        self._resids = np.array(self.block_data.get("auth_seq_id"), dtype=int)
        self._icodes = np.array(self.block_data.get("pdbx_PDB_ins_code"), dtype="<U10")
        _, self._resindices = np.unique(self._resids, return_inverse=True)

    def get_unique_residue_count(self, chains):
        # Combine the chain labels, residue names, and residue ids
        combined = np.core.defchararray.add(chains.auths, self._resnames)  # type: ignore
        combined = np.core.defchararray.add(combined, self._resids.astype(str))  # type: ignore

        return np.unique(combined).size

    @property
    def resnames(self):
        return self._resnames

    @property
    def resids(self):
        return self._resids

    @property
    def resindices(self):
        return self._resindices

    @property
    def icodes(self):
        return self._icodes

    def __len__(self):
        return len(self._resnames)

    def __getitem__(self, index):
        return self._resnames[index], self._resids[index], self._resindices[index]

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


class Atoms:
    def __init__(self, block_data):
        self.block_data = block_data

        self._names = np.array([], dtype="<U10")
        self._ids = np.array([], dtype=int)
        self._elements = np.array([], dtype="<U10")

        self._process()

    def _process(self):
        self._names = np.array(self.block_data.get("label_atom_id"))
        self._ids = np.array(self.block_data.get("id"), dtype=int)
        self._elements = np.core.defchararray.capitalize(
            self.block_data.get("type_symbol")
        )
        self._types = np.array(self.block_data.get("type_symbol"))
        # capitalize the element symbols
        # self._elements = np.core.defchararray.capitalize(self._elements)  # type: ignore

    @property
    def names(self):
        return self._names

    @property
    def types(self):
        return self._types

    @property
    def ids(self):
        return self._ids

    @property
    def elements(self):
        return self._elements

    def __len__(self):
        return len(self._names)

    def __getitem__(self, index):
        return self._names[index], self._ids[index], self._elements[index]

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


class Chains:
    def __init__(self, block_data):
        self.block_data = block_data

        self._labels = np.array([], dtype="<U10")
        self._resindices = np.array([], dtype=int)
        self._ids = np.array([], dtype=int)

        self.mapping = {}

        self._process()

    def _process(self):
        self._labels = np.array(self.block_data.get("label_asym_id"))
        self._auths = np.array(self.block_data.get("auth_asym_id"))
        _, self._ids = np.unique(self._auths, return_inverse=True)
        self._ids += 1

        self.mapping = dict(zip(self._auths, self._ids))

    @property
    def labels(self):
        return self._labels

    @property
    def auths(self):
        return self._auths

    @property
    def ids(self):
        return self._ids

    def __len__(self):
        return len(self._labels)

    def __getitem__(self, index):
        return self._labels[index], self._ids[index]

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]
