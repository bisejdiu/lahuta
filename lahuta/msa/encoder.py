from typing import Tuple

import numpy as np
import pandas as pd
from numpy.typing import NDArray


# pyright: reportUnknownMemberType=false
def encode_labels(
    source_labels: NDArray[np.str_], target_labels: NDArray[np.str_]
) -> Tuple[NDArray[np.int32], NDArray[np.int32]]:
    labels: NDArray[np.str_] = np.concatenate((source_labels.ravel(), target_labels.ravel()))
    indices, _ = pd.factorize(labels)  # type: ignore
    source_pairs = indices[: source_labels.size].reshape(source_labels.shape)
    target_pairs: NDArray[np.int32] = indices[source_labels.size :].reshape(target_labels.shape)

    return source_pairs, target_pairs
