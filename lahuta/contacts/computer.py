import contextlib
import warnings
from typing import Any, Callable, Dict, Generator, Iterable, List, Literal, Optional, Tuple, Union, cast

import joblib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed  # type: ignore
from numpy.typing import NDArray
from tqdm import tqdm
from typing_extensions import TypeAlias

from lahuta.contacts import AtomPlaneContacts, F
from lahuta.core.neighbor_finder import NeighborSearch
from lahuta.core.neighbors import NeighborPairs
from lahuta.core.universe import Universe
from lahuta.lahuta_types.mdanalysis import AtomGroupType

Pairs: TypeAlias = NDArray[np.int32]
Distances: TypeAlias = NDArray[np.float32]
ContactFunction = Callable[[NeighborPairs], NeighborPairs]
ContactFunctions = Union[ContactFunction, Iterable[ContactFunction]]


# pylint: disable=unsubscriptable-object
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
        self, contact_type: Optional[Literal['all', 'atom-atom', 'atom-plane', 'plane-plane']] = 'all'
    ) -> None:
        self._contacts: List[ContactFunction] = []
        if contact_type in ['all', 'atom-atom']:
            for func in self.ATOM_ATOM.values():
                self.register(func)
        if contact_type in ['all', 'atom-plane']:
            self.register_atom_plane_contacts()
        if contact_type in ['all', 'plane-plane']:
            self.register(F.plane_plane_neighbors)

        self.atom_plane_instance: Optional[AtomPlaneContacts] = None
        self._results: Dict[str, NeighborPairs] = {}

    def initialize_atom_plane_instance(self, ns: NeighborPairs) -> None:
        self.atom_plane_instance = AtomPlaneContacts(ns)

    def register_atom_plane_contacts(self):
        self.register([self.donor_pi, self.cation_pi, self.sulphur_pi, self.carbon_pi])

    def donor_pi(self, ns: NeighborPairs) -> NeighborPairs:
        self.initialize_atom_plane_instance(ns)
        assert isinstance(self.atom_plane_instance, AtomPlaneContacts)
        return self.atom_plane_instance.donor_pi()

    def cation_pi(self, ns: NeighborPairs) -> NeighborPairs:
        self.initialize_atom_plane_instance(ns)
        assert isinstance(self.atom_plane_instance, AtomPlaneContacts)
        return self.atom_plane_instance.cation_pi()

    def sulphur_pi(self, ns: NeighborPairs) -> NeighborPairs:
        self.initialize_atom_plane_instance(ns)
        assert isinstance(self.atom_plane_instance, AtomPlaneContacts)
        return self.atom_plane_instance.sulphur_pi()

    def carbon_pi(self, ns: NeighborPairs) -> NeighborPairs:
        self.initialize_atom_plane_instance(ns)
        assert isinstance(self.atom_plane_instance, AtomPlaneContacts)
        return self.atom_plane_instance.carbon_pi()

    def register(self, contact_functions: ContactFunctions) -> None:
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
        if contact_function not in self._contacts:
            raise ValueError(f'{contact_function} is not registered')
        # Unregistering contacts
        self._contacts.remove(contact_function)

    def list_registered(self) -> List[str]:
        # Getting registered contacts
        return [func.__name__ for func in self._contacts]

    def compute(self, ns: NeighborPairs) -> None:
        # Computing all contacts
        self._results = {}
        for func in self._contacts:
            self._results[func.__name__] = func(ns)

    def to_frame(
        self, df_format: Literal["compact", "expanded"] = "expanded", annotations: bool = False
    ) -> pd.DataFrame:
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
        return self._results

    def __getitem__(self, key: str) -> NeighborPairs:
        return self._results[key]

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}(contacts={self.list_registered()})'

    def __str__(self) -> str:
        return self.__repr__()


FrameContacts = Dict[int, Tuple[Pairs, Distances]]


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


@contextlib.contextmanager
def tqdm_joblib(tqdm_object: Any) -> Generator[Any, None, None]:
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""

    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        """
        Class to patch joblib to report into tqdm progress bar given as argument.
        """

        def __call__(self, *args: Any, **kwargs: Any):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)  # type: ignore

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()
