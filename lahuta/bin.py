#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Lahuta CLI
==========
Lahuta is a modern computational tool designed to calculate atom-atom interactions with \
    high performance, scalability, and extensibility.

Authors: Besian I. Sejdiu, PhD
Center for Data Driven Discovery
Department of Structural Biology
Developed at and with ❤️ by St. Jude Children's Research Hospital
"""
import random
from enum import Enum
from time import perf_counter
from typing import Any, Callable, Literal, Optional

import numpy as np
import pandas as pd
from rich import print
from typer import Option, Typer
from typing_extensions import Annotated

from lahuta import Luni
from lahuta.contacts.computer import LahutaContacts
from lahuta.core import NeighborPairs
from lahuta.utils.quote_parser import GROMACS_QUOTES
from lahuta.writers.pretty import PrettyDataFrame

FORMAT = Literal["compact", "expanded"]
RESIDUE_DIFF_DEFAULT = 2
RANDOM_QUOTE = random.choice(GROMACS_QUOTES)


class ContactGroup(str, Enum):
    """Typer Enum group for the supported contact types."""

    all = "all"
    atom_atom = "atom-atom"
    atom_plane = "atom-plane"
    plane_plane = "plane-plane"


class ExportGroup(str, Enum):
    """Typer Enum group for the supported export types."""

    screen = "screen"
    csv = "csv"
    html = "html"
    json = "json"


class FormatGroup(str, Enum):
    """Typer Enum group for the supported dataframe formats."""

    compact = "compact"
    expanded = "expanded"


lahuta = Typer(
    help="[blue]Lahuta CLI.[/blue] :sparkles:",
    rich_markup_mode="rich",
    add_completion=False,
    context_settings={"help_option_names": ["-h", "--help"]},
)


def compute_contacts(ns: NeighborPairs, contact_group: ContactGroup) -> dict[str, NeighborPairs]:
    """Compute contacts.


    Args:
        ns (NeighborPairs): Neighbor pairs.
        contact_group (ContactGroup): Contact group to compute.


    Returns:
        dict[str, NeighborPairs]: Dictionary of contact type and neighbor pairs.

    """
    contacts = LahutaContacts(contact_type=contact_group.value)
    contacts.compute(ns)
    return contacts.results


def compute_neighbors(file_name: str, radius: float, res_dif: int) -> NeighborPairs:
    """Compute neighbors.

    Args:
        file_name (str): File name.
        radius (float): Radius.
        res_dif (int): Residue difference.

    Returns:
        NeighborPairs: Neighbor pairs.

    """
    luni = Luni(file_name)
    return luni.compute_neighbors(radius=radius, res_dif=res_dif)


def create_dataframe(
    result: NeighborPairs, label: str, df_format: FORMAT, n_rows: Optional[int] = None
) -> pd.DataFrame:
    """Create dataframe.

    Args:
        result (NeighborPairs): Neighbor pairs.
        label (str): Label.
        df_format (FORMAT): Dataframe format.
        n_rows (Optional[int], optional): Number of rows. Defaults to None.

    Returns:
        pd.DataFrame: Dataframe.

    """
    if n_rows:
        result = result.clone(result.pairs[:n_rows], result.distances[:n_rows])
    contacts = result.to_frame(df_format=df_format)
    contacts["contact_type"] = np.array([label] * result.pairs.shape[0])
    return contacts


# def dataframe_from_ns(ns: NeighborPairs, df_format: FORMAT) -> pd.DataFrame:


def export_to_screen(results: dict[str, NeighborPairs], df_format: FORMAT, n_rows: int) -> None:
    """Export to screen.

    Args:
        results (dict[str, NeighborPairs]): Dictionary of contact type and neighbor pairs.
        df_format (FORMAT): Dataframe format.
        n_rows (int): Number of rows.

    """
    total_rows = 0
    collect_dfs = []

    for contact_type, result in results.items():
        contact_label = contact_type.split("_")[0]
        row_count = result.pairs.shape[0]
        total_rows += row_count

        if total_rows > n_rows:
            contacts = create_dataframe(result, contact_label, df_format, n_rows)
            collect_dfs.append(contacts)
            break
        else:  # noqa: RET508
            contacts = create_dataframe(result, contact_label, df_format=df_format)
            collect_dfs.append(contacts)

    assert collect_dfs
    total_rows = sum([result.pairs.shape[0] for result in results.values()])
    pretty_df = PrettyDataFrame(pd.concat(collect_dfs), n_rows=n_rows, total_rows=total_rows)
    pretty_df.prettify()


def create_result_dataframe(results: dict[str, NeighborPairs], df_format: FORMAT) -> pd.DataFrame:
    """Create result dataframe.

    Args:
        results (dict[str, NeighborPairs]): Dictionary of contact type and neighbor pairs.
        df_format (FORMAT): Dataframe format.

    Returns:
        pd.DataFrame: Dataframe.

    """
    collect_dfs = []

    for contact_type, result in results.items():
        contact_label = contact_type.split("_")[0]
        contacts = create_dataframe(result, contact_label, df_format)
        collect_dfs.append(contacts)

    return pd.concat(collect_dfs).reset_index(drop=True)


def export_to_csv(results: dict[str, NeighborPairs], df_format: FORMAT) -> None:
    """Export to csv.

    Args:
        results (dict[str, NeighborPairs]): Dictionary of contact type and neighbor pairs.
        df_format (FORMAT): Dataframe format.

    """
    results_df = create_result_dataframe(results, df_format)
    results_df.to_csv("lahuta.csv", index=False)


def export_to_file(
    results: dict[str, NeighborPairs], export_type: str, df_format: FORMAT, output: Optional[str]
) -> None:
    """Export to file.

    Args:
        results (dict[str, NeighborPairs]): Dictionary of contact type and neighbor pairs.
        export_type (str): Export type.
        df_format (FORMAT): Dataframe format.
        output (Optional[str], optional): Output file name. Defaults to None.

    """
    export_functions: dict[str, Callable[[pd.DataFrame, str], None | Any]] = {
        "csv": lambda df, filename: df.to_csv(filename, index=False),
        "html": lambda df, filename: df.to_html(filename),
        "json": lambda df, filename: df.to_json(filename),
    }

    if export_type not in export_functions:
        raise ValueError(f"Unsupported export type: {export_type}")

    results = create_result_dataframe(results, df_format)
    filename = f'{output or "lahuta"}.{export_type}'
    export_functions[export_type](results, filename)


def export_results(
    export: list[ExportGroup],
    results: dict[str, NeighborPairs],
    df_format: FORMAT,
    n_rows: int,
    output: Optional[str] = None,
) -> None:
    """Export results.

    Args:
        export (list[ExportGroup]): list of export types.
        results (dict[str, NeighborPairs]): Dictionary of contact type and neighbor pairs.
        df_format (FORMAT): Dataframe format.
        n_rows (int): Number of rows.
        output (Optional[str], optional): Output file name. Defaults to None.

    """
    for export_type in export:
        if export_type == ExportGroup.screen:
            export_to_screen(results, df_format, n_rows)
        else:
            export_to_file(results, export_type, df_format=df_format, output=output)


@lahuta.command(
    epilog=f'[blue]"{RANDOM_QUOTE[0]}" - {RANDOM_QUOTE[1]} (via GROMACS quotes)[/blue]',
)
def run(
    input: str = Option(..., "--input", "-i", help="Input file name for Luni"),
    radius: Annotated[
        float, Option(..., "--radius", "-r", show_default=True, help="Radius for neighbor computation.")
    ] = 5.0,
    res_dif: Annotated[
        int, Option(..., "--res_dif", "-d", show_default=True, help="Residue difference for neighbor computation.")
    ] = RESIDUE_DIFF_DEFAULT,
    contact_group: ContactGroup = Option(
        "all",
        "--contact-group",
        "-cg",
        help="Contact group to compute.",
    ),
    n_rows: int = Option(20, "--n-rows", "-n", help="Number of rows to display."),
    format: FormatGroup = Option(
        "expanded",
        "--format",
        "-f",
        help="Format of the output.",
    ),
    export: list[ExportGroup] = Option(
        ["screen"],
        "--export",
        "-to",
        help="How to output results.",
    ),
    output: Optional[str] = Option(
        "lahuta",
        "--output",
        "-o",
        help="Name of the output file.",
    ),
) -> None:
    """[blue][bold red]WELCOME TO LAHUTA[/bold red][/blue]

    -----------------

    Authors: Besian I. Sejdiu
    Developed with ❤️ by St. Jude Children's Research Hospital
    Learn more by reading the [link=https://github.com/bisejdiu/lahuta]docs[/link]!

    ---
    Lahuta is a modern computational tool designed to calculate atom-atom interactions with \
high performance, scalability, and extensibility.
    Building on first principles, it provides an intuitive and user-friendly interface that simplifies complex \
computations, delivering reliable results quickly and efficiently.

    While Lahuta is intended to be primarily used through its API, the following \
script/CLI is provided as a convenience for very common operations.

    Should you encounter a bug or want to request a feature, \
please don't hesitate to open an issue at github.com/bisejdiu/lahuta.

    """
    start = perf_counter()
    ns = compute_neighbors(input, radius=radius, res_dif=res_dif)

    results = compute_contacts(ns, contact_group)
    end = perf_counter()

    export_results(export, results, df_format=format.value, n_rows=n_rows, output=output)

    print()
    print(f"Done in: {end - start:.5f} seconds.")


if __name__ == "__main__":
    lahuta()
