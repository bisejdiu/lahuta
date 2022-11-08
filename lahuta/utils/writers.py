"""
Placeholder
"""

from typing import List, Tuple, Union, Dict
from typing_extensions import Protocol, Literal

import numpy as np
import pandas as pd

from ..core.groups import AtomGroup


class DataFrameBase(Protocol):
    """A class for storing and manipulating contact data."""

    def __init__(self, data: pd.DataFrame):
        """A class for storing and manipulating contact data."""

    def dataframe(self) -> pd.DataFrame:
        """Return the dataframe."""


class CompactDataFrame:
    """A class for storing and manipulating contact data."""

    def __init__(self, data: pd.DataFrame):
        """A class for storing and manipulating contact data."""
        self.data = data

    def dataframe(self) -> pd.DataFrame:
        """Build the dataframe based on input format."""

        df = pd.DataFrame.from_dict(self.data)
        df["residue1"] = df.filter(like="residue1_").apply(
            lambda x: "-".join(x.dropna().astype(str)), axis=1
        )
        df["residue2"] = df.filter(like="residue2_").apply(
            lambda x: "-".join(x.dropna().astype(str)), axis=1
        )
        df = df.drop(df.filter(like="residue1_").columns, axis=1)
        df = df.drop(df.filter(like="residue2_").columns, axis=1)

        return df[["residue1", "residue2", "distances"]]

    def __repr__(self):
        return f"<Lahuta DataFrame with {self.data.shape[0]} contacts>"

    def __str__(self):
        return self.__repr__()


class ExpandedDataFrame:
    """A class for storing and manipulating contact data."""

    def __init__(self, data: pd.DataFrame):
        """A class for storing and manipulating contact data."""
        self.data = data

    def dataframe(self) -> pd.DataFrame:
        """Build the dataframe based on input format."""

        return pd.DataFrame.from_dict(self.data)

    def __repr__(self):
        return f"<Lahuta DataFrame with {self.data.shape[0]} contacts>"

    def __str__(self):
        return self.__repr__()


class PrintDataFrame:
    """A class for storing and manipulating contact data."""

    def __init__(self, data: pd.DataFrame):
        """A class for storing and manipulating contact data."""
        self.data = data

    def dataframe(self) -> pd.DataFrame:
        """Build the dataframe based on input format."""

        self.data = {
            k: v for k, v in sorted(self.data.items(), key=lambda item: item[0])
        }
        df = pd.DataFrame.from_dict(self.data)
        df.columns = pd.MultiIndex.from_tuples(
            [tuple(col.split("_")) for col in df.columns]
        )

        return df

    def __repr__(self):
        return f"<Lahuta DataFrame with {self.data.shape[0]} contacts>"

    def __str__(self):
        return self.__repr__()


FACTORY_DICT: Dict[Literal["compact", "expanded", "print"], DataFrameBase] = {
    "compact": CompactDataFrame,
    "expanded": ExpandedDataFrame,
    "print": PrintDataFrame,
}


class DataFrameFactory:
    """A class for storing and manipulating contact data."""

    # TODO:
    # AtomGroup does not neccessarily have a double column format of atoms.
    def __init__(
        self,
        nag: AtomGroup,
        distances: np.ndarray,
        df_format: Literal["compact", "expanded", "print"],
    ):
        """A class for storing and manipulating contact data."""

        col1, col2 = nag[:, 0], nag[:, 1]

        self.methods = ["resids", "resnames", "names", "indices"]
        self.format = df_format

        data = {}
        for method in self.methods:
            data[f"residue1_{method}"] = getattr(col1, method)
            data[f"residue2_{method}"] = getattr(col2, method)
        data["distances"] = distances

        # print(data)

        self.data = pd.DataFrame.from_dict(data)

    def dataframe(self) -> pd.DataFrame:
        """Build the dataframe based on input format."""
        return FACTORY_DICT[self.format](self.data).dataframe()

    def __repr__(self):
        return f"<Lahuta DataFrame with {self.data.shape[0]} contacts>"

    def __str__(self):
        return self.__repr__()
