"""Module to calculate DSSP."""
import hashlib
import io
import shutil
import subprocess
from pathlib import Path
from typing import Any, ClassVar, Optional

import numpy as np
from numpy.typing import NDArray

from lahuta.config._atom_type_strings import BASE_AA_CONVERSION

# class _DSSP:
#     """Class to handle DSSP output.

#     This class provides methods to handle DSSP output.

#     Args:
#         input_file (Optional[str | Path], optional): The input file. Default is None.
#         dssp_dict (Optional[dict[tuple[str, tuple[str, int, str]], tuple[Any, ...]]], optional): The DSSP dictionary.
#         Default is None.

#     Attributes:
#         secondary_structure (NDArray[np.str_]): The secondary structure assignments.
#         solvent_accessible_area (NDArray[np.float32]): The solvent accessible surface areas.

#     Raises:
#         ValueError: If neither input_file nor dssp_dict is provided.
#     """

#     def __init__(
#         self,
#         input_file: Optional[str | Path] = None,
#         dssp_dict: Optional[dict[tuple[str, tuple[str, int, str]], tuple[Any, ...]]] = None,
#     ) -> None:
#         match (input_file, dssp_dict):
#             case (None, None):
#                 raise ValueError("Either input_file or dssp_dict must be provided.")
#             case (None, dssp_dict):
#                 self.dssp_dict = dssp_dict or {}
#             case (input_file, None):
#                 handler = DSSPCLIHandler(str(input_file))
#                 handler.parse_or_calculate_dssp()
#                 self.dssp_dict = handler.dssp_dict

#     @property
#     def secondary_structure(self) -> NDArray[np.str_]:
#         """Return the secondary structure assignments."""
#         # Secondary structure assignments
#         return np.array([data[1] for data in self.dssp_dict.values()])

#     @property
#     def solvent_accessible_area(self) -> NDArray[np.float32]:
#         """Return the solvent accessible surface areas."""
#         # Solvent accessible surface areas
#         return np.array([data[2] for data in self.dssp_dict.values()])


class DSSPParser:
    """Class to parse DSSP output.

    This class can be used to parse DSSP output and return the data as arrays.

    Args:
        dssp_file: The DSSP output file.

    Attributes:
        dssp_file: The DSSP output file.
        aa_conversion: A dictionary to convert amino acid names to one letter codes.
    """

    DTYPES = np.dtype([("chain_auths", "<U2"), ("resname", "<U10"), ("resid", "int")])

    def __init__(self, dssp_file: str | io.StringIO) -> None:
        # dssp_file can be a file path or a file object
        self.dssp_file = dssp_file if hasattr(dssp_file, "read") else Path(dssp_file).open("r")  # noqa: SIM115
        self.aa_conversion = {v: k for k, v in BASE_AA_CONVERSION.items()}

    def parse(
        self,
    ) -> tuple[NDArray[np.void], NDArray[np.str_], NDArray[np.int32]]:
        """Parse the DSSP output and return the data as arrays."""
        resinfo, ss, acc = [], [], []
        start = False
        for line in self.dssp_file:
            if "RESIDUE" in line:
                start = True
                continue
            if not start or len(line) < 136 or line[9] == " ":
                continue

            # Extract relevant data
            res_id = int(line[5:10].strip())
            chain_auth = line[11]
            res_name = line[13]
            sec_struct = line[16].strip() or "-"
            solvent_access = int(line[34:38].strip())

            # Append to lists
            resinfo.append((chain_auth, self.aa_conversion.get(res_name), res_id))
            ss.append(sec_struct)
            acc.append(solvent_access)

        # Convert lists to arrays
        resinfo_array = np.array(resinfo, dtype=DSSPParser.DTYPES)
        ss_array = np.array(ss, dtype="<U1")
        acc_array = np.array(acc, dtype="int")

        return resinfo_array, ss_array, acc_array


class DSSP:
    """Class to handle DSSP command line interface.

    This class can be used to run DSSP from the command line and parse the output.
    The parsed data is stored in the class instance and can be accessed via the
    resinfo_array, ss_array and acc_array attributes.

    Args:
        input_file: Path to the input file.
        executable_name: Name of the DSSP executable. If None, "dssp" is used.

    Attributes:
        input_file: Path to the input file.
        executable_name: Name of the DSSP executable.
        resinfo_array: Array with residue information.
        ss_array: Array with secondary structure information.
        acc_array: Array with solvent accessibility information.
        command: Command to run DSSP.
    """

    _file_hash_cache: ClassVar[dict[str, Any]] = {}

    def __init__(self, input_file: Path | str, executable_name: Optional[str] = None) -> None:
        self.input_file = input_file if isinstance(input_file, Path) else Path(input_file)
        self.executable_name = executable_name
        self._command = self._build_command()
        self._file_hash = self._calculate_file_hash()
        self.resinfo_array: NDArray[np.void] = np.array([], dtype=DSSPParser.DTYPES)
        self.ss_array: NDArray[np.str_] = np.array([], dtype="<U1")
        self.acc_array: NDArray[np.int32] = np.array([], dtype="int")

    def _calculate_file_hash(self) -> str:
        """Calculate a simple hash based on the file name and the first few lines."""
        hasher = hashlib.md5()
        hasher.update(str(self.input_file).encode())

        try:
            with Path.open(self.input_file, "r") as f:
                for _ in range(10):  # Read first 10 lines
                    line = f.readline()
                    hasher.update(line.encode())
        except IOError:
            pass  # If the file cannot be read, use only the name for the hash

        return hasher.hexdigest()

    def parse_or_calculate_dssp(self) -> None:
        """Parse or calculate DSSP output and store the data in arrays."""
        if self._file_hash not in self._file_hash_cache:
            cli_output = self.run()
            file_like_object = io.StringIO(cli_output)

            # Use DSSPParser to parse the output
            parser = DSSPParser(file_like_object)
            self.resinfo_array, self.ss_array, self.acc_array = parser.parse()

            # Cache the result
            self._file_hash_cache[self._file_hash] = (self.resinfo_array, self.ss_array, self.acc_array)
        else:
            # Load from cache
            self.resinfo_array, self.ss_array, self.acc_array = self._file_hash_cache[self._file_hash]

    def _build_command(self) -> list[str]:
        executable_name = self.executable_name if self.executable_name else "dssp"
        return [executable_name, "-i", str(self.input_file)]

    def run(self) -> str:
        """Run DSSP and return the output."""
        if shutil.which(self._command[0]) is None:
            raise RuntimeError(f"Could not find the executable for {self._command[0]}")
        result = subprocess.run(self._command, capture_output=True, text=True, check=True)
        return result.stdout

    @property
    def command(self) -> list[str]:
        """Return the command to run."""
        return self._command
