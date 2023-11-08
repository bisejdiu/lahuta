import io
import shutil
import subprocess
from typing import Any, Optional
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from Bio.PDB.DSSP import _make_dssp_dict


class DSSP:
    def __init__(
        self,
        input_file: Optional[str | Path] = None,
        dssp_dict: Optional[dict[tuple[str, tuple[str, int, str]], tuple[Any, ...]]] = None,
    ) -> None:
        match (input_file, dssp_dict):
            case (None, None):
                raise ValueError('Either input_file or dssp_dict must be provided.')
            case (None, dssp_dict):
                self.dssp_dict = dssp_dict or {}
            case (input_file, None):
                handler = DSSPCLIHandler(str(input_file))
                handler.parse_or_calculate_dssp()
                self.dssp_dict = handler.dssp_dict

    @property
    def secondary_structure(self) -> NDArray[np.str_]:
        # Secondary structure assignments
        return np.array([data[1] for data in self.dssp_dict.values()])

    @property
    def solvent_accessible_area(self) -> NDArray[np.float32]:
        # Solvent accessible surface areas
        return np.array([data[2] for data in self.dssp_dict.values()])


class DSSPCLIHandler:
    def __init__(self, input_file: Path | str, executable_name: Optional[str] = None) -> None:
        self.input_file = input_file
        self.executable_name = executable_name
        self.dssp_dict: dict[tuple[str, tuple[str, int, str]], tuple[Any, ...]] = {}
        self._command = self._build_command()

    def parse_or_calculate_dssp(self) -> None:
        if not self.dssp_dict:
            cli_output = self.run()
            file_like_object = io.StringIO(cli_output)

            dssp_dict, _ = _make_dssp_dict(file_like_object)
            self.dssp_dict = dssp_dict

    def _build_command(self) -> list[str]:
        executable_name = self.executable_name if self.executable_name else 'dssp'
        return [executable_name, '-i', str(self.input_file)]

    def run(self) -> str:
        if shutil.which(self.command[0]) is None:
            raise RuntimeError(f"Could not find the executable for {self.command[0]}")
        result = subprocess.run(self.command, capture_output=True, text=True, check=True)
        return result.stdout

    @property
    def command(self) -> list[str]:
        return self._command


if __name__ == '__main__':
    file_loc = Path(__file__).parent.parent.parent / '1gzm.cif'

    # handler = DSSPCLIHandler(str(file_loc))
    # handler.parse_or_calculate_dssp()
    # dssp = DSSP(dssp_dict=handler.dssp_dict)

    dssp = DSSP(input_file=str(file_loc))

    secondary_structure = dssp.secondary_structure
    print('secondary_structure', secondary_structure)
    # solvent_accessible_area = dssp.solvent_accessible_area
    # print('solvent_accessible_area', solvent_accessible_area)
