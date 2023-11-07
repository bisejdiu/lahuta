import io
import subprocess
from typing import Any

import numpy as np
from numpy.typing import NDArray
from Bio.PDB.DSSP import _make_dssp_dict


class DSSP:
    def __init__(self, dssp_dict: dict[tuple[str, tuple[str, int, str]], tuple[Any, ...]]) -> None:
        self.dssp_dict = dssp_dict

    @property
    def secondary_structure(self) -> NDArray[np.str_]:
        # Secondary structure assignments
        return np.array([data[1] for data in self.dssp_dict.values()])

    @property
    def solvent_accessible_area(self) -> NDArray[np.float32]:
        # Solvent accessible surface areas
        return np.array([data[2] for data in self.dssp_dict.values()])

class DSSPCLIHandler:
    def __init__(self, input_file):
        self.input_file = input_file
        self._command = self._build_command()

    def get_dssp_dict(self):
        cli_output = self.run()
        file_like_object = io.StringIO(cli_output)

        dssp_dict, _ = _make_dssp_dict(file_like_object)
        return dssp_dict

    def _build_command(self):
        return ['mkdssp', '-i', self.input_file]
        
    def run(self):
        result = subprocess.run(self.command, capture_output=True, text=True, check=True)
        return result.stdout

    @property
    def command(self):
        return self._command


if __name__ == '__main__':
    # Test the DSSP class
    # dssp_data = DSSP('/home/bisejdiu/2023/lahuta/1gzm.dssp')
    # secondary_structure = dssp_data.secondary_structure
    # print('secondary_structure', secondary_structure)
    # solvent_accessible_area = dssp_data.solvent_accessible_area
    # print('solvent_accessible_area', solvent_accessible_area)

    handler = DSSPCLIHandler('1gzm.cif')
    dssp_dict = handler.get_dssp_dict()
    dssp_data = DSSP(dssp_dict)
    secondary_structure = dssp_data.secondary_structure
    print('secondary_structure', secondary_structure)
    solvent_accessible_area = dssp_data.solvent_accessible_area
    print('solvent_accessible_area', solvent_accessible_area)
