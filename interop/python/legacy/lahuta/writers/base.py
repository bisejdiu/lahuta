"""Defines the DataFrame abstract base class (ABC).

The DataFrame class establishes an interface for creating pandas DataFrame objects from data 
that is encapsulated in the class. The primary method that subclasses need to implement is 
the `execute` method, which should return a pandas DataFrame object.

"""

from abc import ABC, abstractmethod
from typing import Any

import pandas as pd
from numpy.typing import NDArray


class DataFrame(ABC):
    """Abstract Base Class (ABC) defining the core interface for DataFrame operations.

    This class provides an abstracted layer for the operation of dataframes.
    The 'data' is expected to be a dictionary, where the keys are the column names
    and the values are n-dimensional arrays (ndarrays). The class expects subclasses
    to implement the `execute` method which should return a pandas DataFrame.

    Attributes:
        data (dict[str, NDArray[Any]]): A dictionary mapping column names (str) to
            N-dimensional array-like structures. It represents the data for the dataframe.

    Methods:
        execute: Abstract method which, when implemented by a subclass, should execute
            the necessary dataframe operations and return a pandas DataFrame.
    """

    def __init__(self, data: dict[str, NDArray[Any]]):
        self.data = data

    @abstractmethod
    def execute(self) -> pd.DataFrame:
        """Abstract method that defines the interface for dataframe operations.

        The method is expected to be implemented by any subclass inheriting from this ABC.
        The implementing method should perform necessary dataframe operations and return a
        pandas DataFrame.

        Returns:
            pd.DataFrame: A pandas DataFrame after executing the implemented dataframe operations.

        Raises:
            NotImplementedError: This method should be implemented in a subclass,
                it raises a NotImplementedError if called from an instance of this ABC.
        """
