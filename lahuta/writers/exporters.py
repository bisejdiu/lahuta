"""Module for exporting VMD scripts."""
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

class VMDExporter:
    """Class for exporting VMD scripts.

    Attributes:
        pairs (NDArray[np.int32]): Array of residue pairs.
        TCL_TEMPLATE_PATH (Path): Path to the VMD TCL template file.
        tcl_code (str): VMD TCL code.
    """

    TCL_TEMPLATE_PATH = Path(__file__).parent.parent / "utils" / "vmd_script.tcl"

    def __init__(self, pairs: NDArray[np.int32]):
        self.pairs = pairs
        try:
            with Path.open(self.TCL_TEMPLATE_PATH, "r") as file:
                self.tcl_code = file.read()
        except FileNotFoundError:
            print("COULD NOT LOAD VMD TCL TEMPLATE")

    def export(self, sphere_resolution: int = 20, save_to_file: bool = False) -> str | None:
        """Export VMD script.

        Args:
            sphere_resolution (int, optional): Sphere resolution. Defaults to 20.
            save_to_file (bool, optional): Flag to save the script to a file. Defaults to False.

        Returns:
            str | None: VMD script.
        """
        res_pairs_tcl = "{{{}}}".format(" ".join("{{{0} {1}}}".format(pair[0], pair[1]) for pair in self.pairs))
        tcl_commands = [self.tcl_code, f"draw_interactions top {res_pairs_tcl} {sphere_resolution}"]

        if save_to_file:
            with open("vmd_script.tcl", "w") as file:  # noqa: PTH123
                file.write("\n".join(tcl_commands))
            print("VMD script saved to vmd_script.tcl. Run it: vmd input_file.pdb -e vmd_script.tcl")
            return None

        return "\n".join(tcl_commands)
