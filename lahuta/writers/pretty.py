"""Prettify a pandas dataframe and print it to the console."""
from typing import Any

import pandas as pd
from rich.console import Console
from rich.table import Table


class PrettyDataFrame:
    """Prettify a pandas dataframe and print it to the console.

    Attributes:
        df (pd.DataFrame): The dataframe to be prettified.
        n_rows (int): The number of rows to be printed.
        total_rows (int): The total number of rows in the dataframe.
        console (Console): The rich console object.
        table (Table): The rich table object.

    Methods:
        prettify: Prettify the dataframe and print it to the console.
    """

    INT_COLOR = "red"
    FLOAT_COLOR = "green"
    OTHER_COLOR = "blue"

    def __init__(self, df: pd.DataFrame, n_rows: int = 20, total_rows: int = 0):
        self.total_rows = total_rows
        self.n_rows = n_rows if n_rows < self.total_rows else self.total_rows
        self.df = df.reset_index().rename(columns={"index": ""}).head(self.n_rows)
        self.console = Console()
        self.table = Table(show_footer=False)

    def _add_columns(self) -> None:
        for col in self.df.columns:
            self.table.add_column(str(col))

    def _add_rows(self) -> None:
        for row in self.df.to_numpy():
            formatted_row = self._format_row(row)
            self.table.add_row(*formatted_row)
        if self.n_rows < self.total_rows:
            self.table.add_row(*["..." for _ in range(len(self.df.columns))])

    def _format_row(self, row: list[Any]) -> list[str]:
        formatted_row = []
        for index, item in enumerate(row):
            # Index without color change
            if index == 0:
                formatted_row.append(str(item))
                continue

            if isinstance(item, int):
                formatted_row.append(f"[{self.INT_COLOR}]{item}[/{self.INT_COLOR}]")
            elif isinstance(item, float):
                formatted_row.append(f"[{self.FLOAT_COLOR}]{item:.4f}[/{self.FLOAT_COLOR}]")
            else:
                formatted_row.append(f"[{self.OTHER_COLOR}]{item}[/{self.OTHER_COLOR}]")

        return formatted_row

    def _add_caption(self) -> None:
        self.table.caption = f"Showing {self.n_rows} rows out of {self.total_rows} computed rows"

    def prettify(self) -> None:
        """Prettify the dataframe and print it to the console."""
        self._add_columns()
        self._add_rows()
        self._add_caption()
        self.console.print(self.table)
