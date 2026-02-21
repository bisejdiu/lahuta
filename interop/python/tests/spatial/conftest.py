# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     d = [{"v": "besian"}, {"v": "sejdiu"}, {"v": "@gmail.com"}]
#     print("".join(map(operator.itemgetter("v"), d)))
#
import numpy as np
import pytest


@pytest.fixture
def coords_simple() -> np.ndarray:
    return np.array(
        [
            [0.0, 0.0, 0.0],  # 0
            [1.0, 0.0, 0.0],  # 1
            [0.0, 1.0, 0.0],  # 2
            [0.0, 0.0, 1.0],  # 3
            [1.0, 1.0, 0.0],  # 4
        ],
        dtype=np.float64,
    )
