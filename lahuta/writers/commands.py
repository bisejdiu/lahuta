"""
This module defines two concrete classes, CompactDataFrame and ExpandedDataFrame, 
that inherit from the abstract base class DataFrame.

The CompactDataFrame and ExpandedDataFrame classes provide different implementations 
of the `execute` method inherited from the DataFrame abstract base class.

- CompactDataFrame: The class provides a compact format of a pandas DataFrame. The data 
  is collapsed into fewer columns.
- ExpandedDataFrame: The class provides an expanded format of a pandas DataFrame. The data 
  is returned as is without any manipulation.

These classes are typically used in combination with the DataFrameWriter class, which uses 
a factory design pattern to create either a CompactDataFrame or an ExpandedDataFrame based 
on the provided format string.

Available Classes:
- CompactDataFrame: Concrete class for creating dataframes in compact format.
- ExpandedDataFrame: Concrete class for creating dataframes in expanded format.

Example usage:
    compact_df = CompactDataFrame(data).execute()
    expanded_df = ExpandedDataFrame(data).execute()
"""


from typing import Dict, Type

import pandas as pd

from .base import DataFrame


class CompactDataFrame(DataFrame):
    """Concrete class implementing the DataFrame ABC to create compact dataframes.

    The CompactDataFrame provides an implementation of the `execute` method. It creates a pandas
    DataFrame in a compact format where related columns are collapsed into a single one for each
    partner, and unnecessary columns are dropped.

    Methods:
        execute: Implementation of the `execute` method from the DataFrame ABC. Creates a compact
            format of the pandas DataFrame.
    """

    def execute(self) -> pd.DataFrame:
        """Creates the dataframe in a compact format.

        The method compacts the data by joining the related columns for each partner into a single
        column and dropping unnecessary columns. The resultant dataframe includes "partner1",
        "partner2", and "distances" columns.

        Returns:
            pd.DataFrame: A pandas DataFrame in a compact format.
        """
        df = pd.DataFrame(self.data)
        df["partner1"] = df.filter(like="partner1_").apply(  # type: ignore
            lambda x: "-".join(x.dropna().astype(str)), axis=1
        )
        df["partner2"] = df.filter(like="partner2_").apply(  # type: ignore
            lambda x: "-".join(x.dropna().astype(str)), axis=1
        )
        df = df.drop(df.filter(like="partner1_").columns, axis=1)  # type: ignore
        df = df.drop(df.filter(like="partner2_").columns, axis=1)  # type: ignore

        return df[["partner1", "partner2", "distances"]]


class ExpandedDataFrame(DataFrame):
    """Concrete class implementing the DataFrame ABC for creating expanded dataframes.

    The ExpandedDataFrame provides an implementation of the `execute` method. It simply converts
    the input data into a pandas DataFrame without applying any transformation or manipulation to
    the columns or rows.

    Methods:
        execute: Implementation of the `execute` method from the DataFrame ABC. Returns a pandas
            DataFrame in the same format as the input data.
    """

    def execute(self) -> pd.DataFrame:
        """Creates the dataframe in an expanded format.

        The method simply converts the data into a pandas DataFrame without any manipulation.

        Returns:
            pd.DataFrame: A pandas DataFrame in the same format as the input data.
        """
        return pd.DataFrame(self.data)


FACTORY_DICT: Dict[str, Type[DataFrame]] = {
    "compact": CompactDataFrame,
    "expanded": ExpandedDataFrame,
}
