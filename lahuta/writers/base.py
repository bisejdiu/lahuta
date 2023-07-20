from abc import ABC, abstractmethod

import pandas as pd


class DataFrame(ABC):
    """Abstract base class for dataframe operations."""

    def __init__(self, data):
        self.data = data

    @abstractmethod
    def execute(self) -> pd.DataFrame:
        """Execute the dataframe command."""
