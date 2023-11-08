"""Module to calculate DSSP."""
import io
import shutil
import subprocess
from pathlib import Path
from typing import Any, Optional

import numpy as np
from Bio.PDB.DSSP import _make_dssp_dict
from numpy.typing import NDArray


class DSSP:
    """Class to handle DSSP output.

    This class provides methods to handle DSSP output.

    Args:
        input_file (Optional[str | Path], optional): The input file. Default is None.
        dssp_dict (Optional[dict[tuple[str, tuple[str, int, str]], tuple[Any, ...]]], optional): The DSSP dictionary.
        Default is None.

    Attributes:
        secondary_structure (NDArray[np.str_]): The secondary structure assignments.
        solvent_accessible_area (NDArray[np.float32]): The solvent accessible surface areas.

    Raises:
        ValueError: If neither input_file nor dssp_dict is provided.
    """

    def __init__(
        self,
        input_file: Optional[str | Path] = None,
        dssp_dict: Optional[dict[tuple[str, tuple[str, int, str]], tuple[Any, ...]]] = None,
    ) -> None:
        match (input_file, dssp_dict):
            case (None, None):
                raise ValueError("Either input_file or dssp_dict must be provided.")
            case (None, dssp_dict):
                self.dssp_dict = dssp_dict or {}
            case (input_file, None):
                handler = DSSPCLIHandler(str(input_file))
                handler.parse_or_calculate_dssp()
                self.dssp_dict = handler.dssp_dict

    @property
    def secondary_structure(self) -> NDArray[np.str_]:
        """Return the secondary structure assignments."""
        # Secondary structure assignments
        return np.array([data[1] for data in self.dssp_dict.values()])

    @property
    def solvent_accessible_area(self) -> NDArray[np.float32]:
        """Return the solvent accessible surface areas."""
        # Solvent accessible surface areas
        return np.array([data[2] for data in self.dssp_dict.values()])


class DSSPCLIHandler:
    """Class to handle DSSP command line interface.

    This class provides methods to handle DSSP command line interface.

    Args:
        input_file (Path | str): The input file.
        executable_name (Optional[str], optional): The executable name. Default is None.

    Attributes:
        input_file (Path | str): The input file.
        executable_name (Optional[str]): The executable name.
        dssp_dict (dict[tuple[str, tuple[str, int, str]], tuple[Any, ...]]): The DSSP dictionary.
        command (list[str]): The command to run.

    """

    def __init__(self, input_file: Path | str, executable_name: Optional[str] = None) -> None:
        self.input_file = input_file
        self.executable_name = executable_name
        self.dssp_dict: dict[tuple[str, tuple[str, int, str]], tuple[Any, ...]] = {}
        self._command = self._build_command()

    def parse_or_calculate_dssp(self) -> None:
        """Parse or calculate DSSP output."""
        if not self.dssp_dict:
            cli_output = self.run()
            file_like_object = io.StringIO(cli_output)

            dssp_dict, _ = _make_dssp_dict(file_like_object)
            self.dssp_dict = dssp_dict

    def _build_command(self) -> list[str]:
        executable_name = self.executable_name if self.executable_name else "dssp"
        return [executable_name, "-i", str(self.input_file)]

    def run(self) -> str:
        """Run DSSP."""
        if shutil.which(self.command[0]) is None:
            raise RuntimeError(f"Could not find the executable for {self.command[0]}")
        result = subprocess.run(self.command, capture_output=True, text=True, check=True)
        return result.stdout

    @property
    def command(self) -> list[str]:
        """Return the command to run."""
        return self._command
