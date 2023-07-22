from abc import ABC, abstractmethod
from typing import Any, Dict

import pandas as pd
from numpy.typing import NDArray


class DataFrame(ABC):
    """Abstract base class for dataframe operations."""

    def __init__(self, data: Dict[str, NDArray[Any]]):
        self.data = data

    @abstractmethod
    def execute(self) -> pd.DataFrame:
        """Execute the dataframe command."""
