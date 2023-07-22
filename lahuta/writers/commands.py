from typing import Dict, Type

import pandas as pd

from .base import DataFrame


class CompactDataFrame(DataFrame):
    """Concrete class for creating compact dataframes."""

    def execute(self) -> pd.DataFrame:
        """Create the dataframe in compact format."""
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
    """Concrete class for creating expanded dataframes."""

    def execute(self) -> pd.DataFrame:
        """Create the dataframe in expanded format."""
        return pd.DataFrame(self.data)


FACTORY_DICT: Dict[str, Type[DataFrame]] = {
    "compact": CompactDataFrame,
    "expanded": ExpandedDataFrame,
}
