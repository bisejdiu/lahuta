from typing import Any

import numpy as np
from numpy.typing import NDArray
from Bio.PDB.DSSP import make_dssp_dict


class DSSPData:
    def __init__(self, dssp_file_path: str) -> None:
        self.dssp_file_path: str = dssp_file_path
        self.dssp_dict: dict[tuple[str, tuple[str, int, str]], tuple[Any, ...]] = self._parse_dssp_file()

    def _parse_dssp_file(self) -> dict[tuple[str, tuple[str, int, str]], tuple[Any, ...]]:
        # Call the make_dssp_dict function from BioPython and return only the dssp_dict
        dssp_dict, _ = make_dssp_dict(self.dssp_file_path)
        return dssp_dict  # type: ignore

    @property
    def secondary_structure(self) -> NDArray[np.str_]:
        # Secondary structure assignments
        return np.array([data[1] for data in self.dssp_dict.values()])

    @property
    def solvent_accessible_area(self) -> NDArray[np.float32]:
        # Solvent accessible surface areas
        return np.array([data[2] for data in self.dssp_dict.values()])


if __name__ == '__main__':
    # Test the DSSPData class
    dssp_data = DSSPData('/home/bisejdiu/2023/lahuta/1gzm.dssp')
    secondary_structure = dssp_data.secondary_structure
    print('secondary_structure', secondary_structure)
    solvent_accessible_area = dssp_data.solvent_accessible_area
    print('solvent_accessible_area', solvent_accessible_area)
