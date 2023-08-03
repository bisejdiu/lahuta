import warnings
from typing import Callable, Dict, Iterable, List, Literal, Optional, Tuple, Union, cast

import numpy as np
import pandas as pd
from joblib import Parallel, delayed  # type: ignore
from MDAnalysis.analysis.base import AnalysisBase
from numpy.typing import NDArray
from tqdm import tqdm
from typing_extensions import TypeAlias

from lahuta.contacts import AtomPlaneContacts, F
from lahuta.core.neighbor_finder import NeighborSearch
from lahuta.core.neighbors import NeighborPairs
from lahuta.core.universe import Universe
from lahuta.lahuta_types.mdanalysis import AtomGroupType

from ._ctx_mngrs import tqdm_joblib

Pairs: TypeAlias = NDArray[np.int32]
Distances: TypeAlias = NDArray[np.float32]
ContactFunction = Callable[[NeighborPairs], NeighborPairs]
ContactFunctions = Union[ContactFunction, Iterable[ContactFunction]]
FrameContacts = Dict[int, Tuple[Pairs, Distances]]


# pylint: disable=unsubscriptable-object
class LahutaContacts:
    """
    A class to manage the computation of various molecular contacts.

    LahutaContacts is a class designed to manage the computation of various molecular contacts.
    It allows users to specify the types of contacts they want to compute, such as atom-atom,
    atom-plane, or plane-plane contacts, and provides methods to register, unregister, and
    compute those contacts.

    Methods:
        register(contact_functions): Registers given contact functions.
        unregister(contact_function): Unregisters the given contact function.
        list_registered(): Returns the list of registered contact function names.
        compute(ns): Computes all the registered contacts.
        to_frame(df_format="expanded", annotations=False): Converts the computed contacts to a DataFrame.
        results: Property to access the computed results.

    Examples:
        ns = NeighborPairs(...)

        contacts = LahutaContacts(contact_type='atom-atom') # or 'all', 'atom-plane', 'plane-plane' or None
        contacts.register(F.custom_contact_function) # custom_contact_function must be a callable function
        contacts.unregister(F.vdw_neighbors) # example of unregistering a function (function must be registered first)
        contacts.list_registered() # returns a list of registered functions
        contacts.compute(ns)
        results = contacts.results
        result_df = contacts.to_frame(annotations=True)

    Note:
        - Users can register custom-defined contact functions.
        - Default contact_type is 'all', which will include all predefined contact types.
        - Warnings are generated if attempts are made to register functions that are already registered, or
            if non-callable objects are attempted to be registered.
        - A ValueError is raised if an attempt is made to unregister a function that is not already registered.
    """

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
        self, contact_type: Optional[Literal['all', 'atom-atom', 'atom-plane', 'plane-plane']] = 'all'
    ) -> None:
        self._contacts: List[ContactFunction] = []
        if contact_type in ['all', 'atom-atom']:
            for func in self.ATOM_ATOM.values():
                self.register(func)
        if contact_type in ['all', 'atom-plane']:
            self._register_atom_plane_contacts()
        if contact_type in ['all', 'plane-plane']:
            self.register(F.plane_plane_neighbors)

        self.atom_plane_instance: Optional[AtomPlaneContacts] = None
        self._results: Dict[str, NeighborPairs] = {}

    def _initialize_atom_plane_instance(self, ns: NeighborPairs) -> None:
        self.atom_plane_instance = AtomPlaneContacts(ns)

    def _register_atom_plane_contacts(self):
        self.register([self.donor_pi, self.cation_pi, self.sulphur_pi, self.carbon_pi])

    def donor_pi(self, ns: NeighborPairs) -> NeighborPairs:
        """Wrapper for AtomPlaneContacts.donor_pi"""
        self._initialize_atom_plane_instance(ns)
        assert isinstance(self.atom_plane_instance, AtomPlaneContacts)
        return self.atom_plane_instance.donor_pi()

    def cation_pi(self, ns: NeighborPairs) -> NeighborPairs:
        """Wrapper for AtomPlaneContacts.cation_pi"""
        self._initialize_atom_plane_instance(ns)
        assert isinstance(self.atom_plane_instance, AtomPlaneContacts)
        return self.atom_plane_instance.cation_pi()

    def sulphur_pi(self, ns: NeighborPairs) -> NeighborPairs:
        """Wrapper for AtomPlaneContacts.sulphur_pi"""
        self._initialize_atom_plane_instance(ns)
        assert isinstance(self.atom_plane_instance, AtomPlaneContacts)
        return self.atom_plane_instance.sulphur_pi()

    def carbon_pi(self, ns: NeighborPairs) -> NeighborPairs:
        """Wrapper for AtomPlaneContacts.carbon_pi"""
        self._initialize_atom_plane_instance(ns)
        assert isinstance(self.atom_plane_instance, AtomPlaneContacts)
        return self.atom_plane_instance.carbon_pi()

    def register(self, contact_functions: ContactFunctions) -> None:
        """
        Register contact functions.

        This method takes in a single contact function or an iterable of contact functions and registers them.

        Args:
            contact_functions (ContactFunctions): A single contact function or an iterable of contact functions.

        Returns:
            None

        Raises:
            ValueError: If a non-callable object is passed in.

        """
        # Ensure contact_functions is an iterable
        if not isinstance(contact_functions, Iterable):
            contact_functions = [contact_functions]

        for func in contact_functions:
            if callable(func):
                if func in self._contacts:
                    warnings.warn(f"{func.__name__} is already registered. Skipping.")
                else:
                    self._contacts.append(func)
            else:
                warnings.warn(f"{func} is not a callable function and cannot be registered.")

    def unregister(self, contact_function: ContactFunction) -> None:
        """
        Unregister contact functions.

        This method takes in a single contact function and unregisters it.

        Args:
            contact_function (ContactFunction): A single contact function.

        Returns:
            None

        Raises:
            ValueError: If the contact function is not registered.

        """
        if contact_function not in self._contacts:
            raise ValueError(f'{contact_function} is not registered')
        # Unregistering contacts
        self._contacts.remove(contact_function)

    def list_registered(self) -> List[str]:
        """
        List registered contact functions.

        Args:
            None

        Returns:
            List[str]: A list of the names of the registered contact functions.

        """

        # Getting registered contacts
        return [func.__name__ for func in self._contacts]

    def compute(self, ns: NeighborPairs) -> None:
        """
        Compute all registered contacts.

        This method takes in a NeighborPairs object and computes all the registered contacts.

        Args:
            ns (NeighborPairs): A NeighborPairs object.

        Returns:
            None

        """
        # Computing all contacts
        self._results = {}
        for func in self._contacts:
            self._results[func.__name__] = func(ns)

    def to_frame(
        self, df_format: Literal["compact", "expanded"] = "expanded", annotations: bool = False
    ) -> pd.DataFrame:
        """
        Convert computed contacts to a DataFrame.

        Args:
            df_format (Literal["compact", "expanded"], optional): The format of the DataFrame. Default is "expanded".
            annotations (bool, optional): Whether to include annotations. Default is False.

        Returns:
            pd.DataFrame: A DataFrame containing the computed contacts.

        """
        df = pd.DataFrame()  # type: ignore
        for key, value in self._results.items():
            df_contact = value.to_frame(df_format=df_format, annotations=annotations)
            df_contact['contact_type'] = key
            if key == 'plane_plane_neighbors':
                df_contact_type = value.to_frame(df_format='expanded', annotations=True)['contact_labels']  # type: ignore
                df_contact['contact_labels'] = df_contact_type

            df = pd.concat([df, df_contact], axis=0)  # type: ignore
        return df

    @property
    def results(self) -> Dict[str, NeighborPairs]:
        """
        Return the computed results.
        """
        return self._results

    def __getitem__(self, key: str) -> NeighborPairs:
        return self._results[key]

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}(contacts={self.list_registered()})'

    def __str__(self) -> str:
        return self.__repr__()


class LahutaTrajectoryContacts:
    def __init__(self, res_dif: int, radius: float):
        self.res_dif = res_dif
        self.radius = radius
        self.results: Dict[int, Union[NeighborPairs, Dict[str, NeighborPairs]]] = {}

    def _get_block_slices(self, n_frames: int, n_blocks: int):
        n_frames_per_block = n_frames // n_blocks
        blocks = [range(i * n_frames_per_block, (i + 1) * n_frames_per_block) for i in range(n_blocks - 1)]
        blocks.append(range((n_blocks - 1) * n_frames_per_block, n_frames))
        return blocks

    def _per_frame_compute(self, mda: AtomGroupType, frame_index: int, result: FrameContacts) -> FrameContacts:
        neighbors = NeighborSearch(mda)
        pairs, distances = neighbors.compute(
            radius=self.radius,
            res_dif=self.res_dif,
        )

        result[frame_index] = (pairs, distances)
        return result

    def _frame_iter(self, mda: AtomGroupType, blockslice: slice) -> FrameContacts:
        result: FrameContacts = {}
        trajectory = mda.universe.trajectory
        for ts in trajectory[blockslice.start : blockslice.stop]:
            result = self._per_frame_compute(mda, ts.frame, result)

        return result

    def _compute(self, mda: AtomGroupType, n_jobs: int) -> List[FrameContacts]:
        n_frames = mda.universe.trajectory.n_frames
        blocks = self._get_block_slices(n_frames, n_blocks=n_jobs)
        with tqdm_joblib(tqdm(total=n_jobs)) as _:
            results: List[FrameContacts] = cast(
                List[FrameContacts], Parallel(n_jobs=n_jobs)(delayed(self._frame_iter)(mda, bs) for bs in blocks)
            )

        return results

    def compute(self, luni: Universe, lahuta_contacts: Optional["LahutaContacts"] = None, n_jobs: int = 1):
        if not self.results:
            self.results = {}

        mda = luni.to("mda")
        # get the very first result from results
        results = self._compute(mda, n_jobs=n_jobs)

        ref_ns = None
        for block in results:
            for frame_index, (pairs, distances) in block.items():
                # Check if this is the first frame to create ref_ns
                if ref_ns is None:
                    ref_ns = NeighborPairs(mda, luni.to("mol"), luni.atom_types, pairs, distances)

                ns = ref_ns.clone(pairs, distances)
                if lahuta_contacts is None:
                    self.results[frame_index] = ns
                else:
                    lahuta_contacts.compute(ns)
                    self.results[frame_index] = lahuta_contacts.results


class SlowLahutaTrajectoryContacts(AnalysisBase):
    def __init__(self, luni: Universe, res_dif: int, radius: float):
        self.luni = luni
        self.res_dif = res_dif
        self.radius = radius
        self._trajectory = self.luni.to("mda").universe.trajectory
        self.results: Dict[int, Union[NeighborPairs, Dict[str, NeighborPairs]]] = {}
        self.lahuta_contacts: Optional["LahutaContacts"] = None
        AnalysisBase.__init__(self, self._trajectory)  # type: ignore

    def _single_frame(self):
        """
        Compute the contacts for a single frame.
        """
        ns = self.luni.compute_neighbors(res_dif=self.res_dif, radius=self.radius)
        if self.lahuta_contacts is None:
            self.results[self._frame_index] = ns
        else:
            self.lahuta_contacts.compute(ns)
            self.results[self._frame_index] = self.lahuta_contacts.results

    def compute(self, lahuta_contacts: Optional["LahutaContacts"] = None):
        """
        Compute the contacts for the whole trajectory.
        """
        self.lahuta_contacts = lahuta_contacts
        self.run()  # type: ignore
