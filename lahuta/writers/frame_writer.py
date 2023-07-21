from typing import TYPE_CHECKING, Any, Dict, Literal, Optional, Sequence

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from lahuta.core.neighbors import NeighborPairs

from .commands import FACTORY_DICT


class DataFrameWriter:
    """Class for creating dataframe operations."""

    def __init__(
        self,
        ns: "NeighborPairs",
        df_format: Literal["compact", "expanded"] = "expanded",
        annotations: Optional[Dict[str, Sequence[Any]]] = None,
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
        """Create the dataframe operation."""
        data = self.build()
        if self.format == "compact":
            return FACTORY_DICT[self.format](data).execute()
        elif self.format == "expanded":
            return FACTORY_DICT[self.format](data).execute()
        else:
            raise ValueError("Format must be compact or expanded.")

    def build(self):
        """Build the data."""
        attrs = ["resids", "resnames", "names", "indices"]
        p1, p2 = self.ns.partner1, self.ns.partner2
        data = {f"partner1_{method}": getattr(p1, method) for method in attrs}
        data.update({f"partner2_{method}": getattr(p2, method) for method in attrs})
        data["distances"] = self.ns.distances

        # if annotations is not None, add them to data
        if self.annotations:
            data.update(self.annotations)

        return data
