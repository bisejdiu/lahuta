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
