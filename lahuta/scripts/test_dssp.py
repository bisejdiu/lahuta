from pathlib import Path
from lahuta.analysis.dssp import DSSP
from lahuta.tests.data import Rhodopsin

if __name__ == '__main__':
    file_loc = Rhodopsin().file_loc
    # file_loc = Path(__file__).parent.parent.parent / '1gzm.cif'

    # handler = DSSPCLIHandler(str(file_loc))
    # handler.parse_or_calculate_dssp()
    # dssp = DSSP(dssp_dict=handler.dssp_dict)

    dssp = DSSP(input_file=file_loc)
    dssp.parse_or_calculate_dssp()

    # secondary_structure = dssp.secondary_structure
    secondary_structure = dssp.ss_array
    print('secondary_structure', secondary_structure.shape)
    # solvent_accessible_area = dssp.solvent_accessible_area
    # print('solvent_accessible_area', solvent_accessible_area)
