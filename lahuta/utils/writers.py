"""
Placeholder
"""

from typing import Dict, Any, Protocol, Literal

import numpy as np
import pandas as pd

from ..core.groups import AtomGroup

from dataclasses import dataclass


@dataclass
class DataFrameBase(Protocol):
    """A class for storing and manipulating contact data."""

    def dataframe(self, data: Dict[Any, Any]) -> pd.DataFrame:
        ...


@dataclass
class CompactDataFrame:
    """A class for storing and manipulating contact data."""

    def dataframe(self, data: Dict[Any, Any]) -> pd.DataFrame:
        """Build the dataframe based on input format."""

        self.data = data
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


@dataclass
class PrintDataFrame:
    """A class for storing and manipulating contact data."""

    def dataframe(self, data: Dict[Any, Any]) -> pd.DataFrame:
        """Build the dataframe based on input format."""

        data = {k: v for k, v in sorted(data.items(), key=lambda item: item[0])}
        self.data = pd.DataFrame.from_dict(data)
        self.data.columns = pd.MultiIndex.from_tuples(
            [tuple(col.split("_")) for col in self.data.columns]
        )

        return self.data


@dataclass
class ExpandedDataFrame:
    """A class for storing and manipulating contact data."""

    def dataframe(self, data: Dict[Any, Any]) -> pd.DataFrame:
        """Build the dataframe based on input format."""
        self.data = data

        return pd.DataFrame.from_dict(data)


class DataFrameFactory:
    """A class for storing and manipulating contact data."""

    # TODO:
    # AtomGroup does not neccessarily have a double column format of atoms.
    def __init__(
        self,
        nag: AtomGroup,
        distances: np.ndarray,
        df_format: Literal["print", "compact", "expanded"] = "print",
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

        self.data = data

    def dataframe(self) -> pd.DataFrame:
        """Build the dataframe based on input format."""

        return FACTORY_DICT[self.format].dataframe(self.data)


FACTORY_DICT: Dict[str, DataFrameBase] = {
    "compact": CompactDataFrame(),
    "expanded": ExpandedDataFrame(),
    "print": PrintDataFrame(),
}
