from typing import Callable, Dict, Iterable, List, Literal, Optional, Union

import pandas as pd

from lahuta.contacts import AtomPlaneContacts, F
from lahuta.core.neighbors import NeighborPairs

ContactFunction = Callable[[NeighborPairs], NeighborPairs]
ContactFunctions = Union[ContactFunction, Iterable[ContactFunction]]


class LahutaContacts:
    ATOM_ATOM: Dict[str, ContactFunction] = {
        'aromatic': F.aromatic_neighbors,
        'ionic': F.ionic_neighbors,
        'carbonyl': F.carbonyl_neighbors,
        'vdw': F.vdw_neighbors,
        'hydrophobic': F.hydrophobic_neighbors,
        'hbond': F.hbond_neighbors,
        'weak_hbond': F.weak_hbond_neighbors,
        'polar_hbond': F.polar_hbond_neighbors,
        'weak_polar_hbond': F.weak_polar_hbond_neighbors,
        'metalic': F.metalic_neighbors,
    }

    def __init__(
        self,
        ns: NeighborPairs,
        contact_type: Optional[Literal['all', 'atom-atom', 'atom-plane', 'plane-plane']] = 'all',
    ):
        self.ns = ns
        self._contacts: List[Callable[[NeighborPairs], NeighborPairs]] = []
        if contact_type in ['all', 'atom-atom']:
            for func in self.ATOM_ATOM.values():
                self.register(func)
        if contact_type in ['all', 'atom-plane']:
            self.register_atom_plane_contacts()
        if contact_type in ['all', 'plane-plane']:
            self.register(F.plane_plane_neighbors)

        self._results: Dict[str, NeighborPairs] = {}

    def register_atom_plane_contacts(self):
        # Wrapper method to register atom-plane contacts
        self.atom_plane_instance = AtomPlaneContacts(self.ns)
        self.register([self._donor_pi, self._cation_pi, self._sulphur_pi, self._carbon_pi])

    def _donor_pi(self, _: NeighborPairs) -> NeighborPairs:
        return self.atom_plane_instance.donor_pi()

    def _cation_pi(self, _: NeighborPairs) -> NeighborPairs:
        return self.atom_plane_instance.cation_pi()

    def _sulphur_pi(self, _: NeighborPairs) -> NeighborPairs:
        return self.atom_plane_instance.sulphur_pi()

    def _carbon_pi(self, _: NeighborPairs) -> NeighborPairs:
        return self.atom_plane_instance.carbon_pi()

    def register(self, contact_functions: ContactFunctions):
        # Registering contacts
        if isinstance(contact_functions, Iterable):
            self._contacts.extend(contact_functions)
        elif callable(contact_functions):
            self._contacts.append(contact_functions)

    def unregister(self, contact_functions: Callable[[NeighborPairs], NeighborPairs]):
        self._contacts.remove(contact_functions)

    def list_registered(self):
        return [func.__name__.lstrip('_') for func in self._contacts]

    def compute(self):
        # Computing all contacts
        self._results = {func.__name__.lstrip('_'): func(self.ns) for func in self._contacts}

    def to_frame(self, df_format: Literal["compact", "expanded"] = "expanded", annotations: bool = False):
        df = pd.DataFrame()
        for key, value in self._results.items():
            df_contact = value.to_frame(df_format=df_format, annotations=annotations)
            df_contact['contact_type'] = key.lstrip('_')
            if key == 'plane_plane_neighbors':
                df_contact_type = value.to_frame(df_format='expanded', annotations=True)['contact_labels']  # type: ignore
                df_contact['contact_labels'] = df_contact_type

            df = pd.concat([df, df_contact], axis=0)  # type: ignore
        return df

    @property
    def results(self):
        return self._results

    def __getitem__(self, key: str):
        return self._results[key]

    def __repr__(self):
        return f'{self.__class__.__name__}({self.ns})'

    def __str__(self):
        return self.__repr__()
