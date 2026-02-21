# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         def __init__(self):
#             self.__class__.r = "besian" + "sejdiu" + "@gmail.com"
#     print(Email().r)
#
from __future__ import annotations

from typing import Iterable

import numpy as np

from .lib.lahuta import KDIndex  # low-level

ArrayLike = np.ndarray | Iterable[Iterable[float]]


# fmt: off
class NearestNeighbors:
    """Minimal sklearn-like wrapper around Lahuta neighbors.

    Parameters:
        radius: Search radius.
        algorithm: 'kd_tree' (default) or 'grid'.
            - 'kd_tree':  builds a KD index on `fit(X)` and uses it for cross queries.
            - 'grid':     use compiled functional path for self and cross. No persistent index.
        sort_results: Whether to sort neighbors by distance per query (costs extra).

    Behavior and performance characteristics:
        - Self queries (`X is None`) always use the optimized compiled self-search path
          (grid based) and exclude the identity element, matching scikit-learn's behavior.
        - Cross queries use the KD index if available (`algorithm='kd_tree'` and `fit` was
          called), otherwise fall back to the compiled functional path which internally builds a
          KD index per call.
    """

    def __init__(self, radius: float, algorithm: str = "kd_tree", sort_results: bool = False) -> None:
        self.radius = float(radius)
        self.algorithm = algorithm
        self.sort_results = bool(sort_results)
        self._X: np.ndarray | None = None
        self._kd: KDIndex | None = None

    def fit(self, X: ArrayLike) -> "NearestNeighbors":
        X = np.asarray(X)
        if X.ndim != 2 or X.shape[1] != 3:
            raise ValueError("X must have shape (n, 3)")
        self._X = X
        if self.algorithm.startswith("kd"):
            kd = KDIndex()
            # Fast paths for zero copy KD build
            if X.flags.c_contiguous and X.dtype in (np.float32, np.float64):
                kd.build_view(X)
            else:
                # Fall back: copy to float64 C and build. There is no dtype conversion in KD
                kd.build(np.asarray(X, dtype=np.float64, order="C"))
            self._kd = kd
        else:
            self._kd = None
        return self

    def radius_neighbors(
        self,
        X: ArrayLike | None = None,
        *,
        return_distance: bool = False,
        sort_results: bool | None = None,
    ) -> list[np.ndarray] | tuple[list[np.ndarray], list[np.ndarray]]:
        if self._X is None:
            raise RuntimeError("Call fit(X) before radius_neighbors().")
        sort = self.sort_results if sort_results is None else bool(sort_results)

        if X is None:
            # Self neighbors: always use grid-based self search
            from .lib.lahuta import neighbors as _cpp_neighbors

            return _cpp_neighbors.radius_neighbors(
                self._X, self.radius, return_distance=return_distance, sort_results=sort
            )

        Q = np.asarray(X, dtype=np.float64)
        if Q.ndim != 2 or Q.shape[1] != 3:
            raise ValueError("X must have shape (m, 3)")

        if self._kd is None or not self.algorithm.startswith("kd"):
            # Functional cross API uses KD internally, but rebuilds each call (no persistent index)
            from .lib.lahuta import neighbors as _cpp_neighbors

            return _cpp_neighbors.radius_neighbors(
                Q, self.radius, Y=self._X, return_distance=return_distance, sort_results=sort
            )

        return self._kd.radius_search(
            Q,
            self.radius,
            grouped=True,
            return_distance=return_distance,
            sort_results=sort,
        )


__all__ = ["NearestNeighbors"]
