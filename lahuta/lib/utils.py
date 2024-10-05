from collections import namedtuple
from typing import overload

import numpy as np
from numpy.typing import NDArray

from lahuta.lib import cLuni, factorize_residues

Labels = list[str] | NDArray[np.str_]
Indices = list[int] | NDArray[np.int32]

FactorizedResidues = namedtuple("FactorizedResidues", ["resindices", "resnames", "resids", "chains"])


@overload
def factorize(resnames: Labels) -> Indices: ...


@overload
def factorize(resnames: Labels, resids: Indices, chainlabels: Labels) -> FactorizedResidues: ...


def factorize(
    resnames: Labels, resids: None | Indices = None, chainlabels: None | Labels = None
) -> Indices | FactorizedResidues:
    if resids is None and chainlabels is None:
        return cLuni.factorize(resnames)
    if resids is not None and chainlabels is not None:
        result = factorize_residues(resnames, resids, chainlabels)
        return FactorizedResidues(*result)

    raise TypeError(
        "Invalid arguments provided. Provide either just resnames or all three: resnames, resids, and chainlabels."
    )
