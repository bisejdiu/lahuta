"""
This module defines the DataFrameWriter class, which serves as a factory for creating 
dataframes in either 'compact' or 'expanded' format from instances of the `NeighborPairs` 
class. 

The DataFrameWriter class is designed to generate dataframes in a specific format, based 
on the provided 'format' argument, from a `NeighborPairs` instance. The class also supports 
the addition of extra annotations to the dataframe. 

Available Classes:
- DataFrameWriter: Factory class for creating DataFrame operations.

Example usage:
    df_writer = DataFrameWriter(ns_instance, df_format='compact')
    df = df_writer.create()
"""

from typing import TYPE_CHECKING, Any, Dict, Literal, Optional

import pandas as pd
from numpy.typing import NDArray

from lahuta.writers.commands import FACTORY_DICT

if TYPE_CHECKING:
    from lahuta.core.neighbors import NeighborPairs


class DataFrameWriter:
    """Class for creating DataFrame operations in either compact or expanded formats.

    This class serves as a factory for creating dataframes in specific formats from
    `NeighborPairs` instances. The formats can be either 'compact' or 'expanded', which
    correspond to CompactDataFrame and ExpandedDataFrame respectively. Additional annotations
    can be added as well, which will be included in the final dataframe.

    Attributes:
        ns (NeighborPairs): Instance of NeighborPairs containing partner information to be
            transformed into a DataFrame.
        format (Literal["compact", "expanded"]): The format in which to create the DataFrame.
            Must be either 'compact' or 'expanded'. Defaults to 'expanded'.
        annotations (Optional[Dict[str, NDArray[Any]]]): Optional additional annotations to
            be included in the DataFrame. The annotations should be a dictionary mapping column
            names (str) to N-dimensional array-like structures.

    Methods:
        create: Constructs the DataFrame in the requested format.
        build: Assembles the data dictionary to be used for DataFrame construction.
    """

    def __init__(
        self,
        ns: "NeighborPairs",
        df_format: Literal["compact", "expanded"] = "expanded",
        annotations: Optional[Dict[str, NDArray[Any]]] = None,
    ):
        """Initialize the factory with a builder and a format."""
        self.ns = ns
        # assert df_format in [
        #     "compact",
        #     "expanded",
        # ], "Format must be compact or expanded."
        self.format = df_format
        self.annotations = annotations

    def create(self) -> pd.DataFrame:
        """Creates the DataFrame operation in the requested format.

        The method first builds the data using the `build` method and then creates a DataFrame
        in the requested format. If the format is not 'compact' or 'expanded', a ValueError is raised.

        Returns:
            pd.DataFrame: The constructed DataFrame in the requested format.

        Raises:
            ValueError: If the format is not 'compact' or 'expanded'.
        """
        data = self.build()
        if self.format == "compact":
            return FACTORY_DICT[self.format](data).execute()
        elif self.format == "expanded":
            return FACTORY_DICT[self.format](data).execute()
        else:
            raise ValueError("Format must be compact or expanded.")

    def build(self) -> Dict[str, NDArray[Any]]:
        """Assembles the data to be used for DataFrame construction.

        The method constructs a dictionary where the keys are column names and the values are the
        corresponding data arrays. It includes partner information and distances from the NeighborPairs
        instance and optional additional annotations.

        Returns:
            Dict[str, NDArray[Any]]: The assembled data to be transformed into a DataFrame.
        """
        attrs = ["resids", "resnames", "names", "indices"]
        p1, p2 = self.ns.partner1, self.ns.partner2
        data = {f"partner1_{method}": getattr(p1, method) for method in attrs}
        data.update({f"partner2_{method}": getattr(p2, method) for method in attrs})
        data["distances"] = self.ns.distances

        # if annotations is not None, add them to data
        if self.annotations:
            data.update(self.annotations)

        return data
