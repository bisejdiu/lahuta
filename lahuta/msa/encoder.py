"""Encoder for MSA labels."""

import numpy as np
import pandas as pd
from numpy.typing import NDArray


def encode_labels(
    source_labels: NDArray[np.void], target_labels: NDArray[np.void]
) -> tuple[NDArray[np.int32], NDArray[np.int32]]:
    """Encode the labels of the source and target sequences.

    Args:
        source_labels (NDArray[np.void]): The labels of the source sequence.
        target_labels (NDArray[np.void]): The labels of the target sequence.

    Returns:
        tuple[NDArray[np.int32], NDArray[np.int32]]: The encoded labels of the source and target sequences.

    """
    labels: NDArray[np.str_] = np.concatenate((source_labels.ravel(), target_labels.ravel()))
    indices, _ = pd.factorize(labels)
    source_pairs = indices[: source_labels.size].reshape(source_labels.shape)
    target_pairs: NDArray[np.int32] = indices[source_labels.size :].reshape(target_labels.shape)

    return source_pairs, target_pairs
