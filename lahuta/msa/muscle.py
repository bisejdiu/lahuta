"""MUSCLE wrapper."""
import subprocess
from typing import Optional


class MuscleFailed(Exception):
    """Exception raised when MUSCLE fails to align the sequences."""


class MuscleCommand:
    """Helper class to construct and execute the MUSCLE command.

    Attributes:
        input_file (str): Path to the input file.
        output_file (str): Path to the output file.
        strategy (Optional[str]): Advanced alignment strategy.
        replicates (Optional[int]): Number of replicates.
        perm (Optional[str]): Permutation method.
        perturb (Optional[int]): Number of perturbations.
        results (Optional[str]): Results of the MUSCLE command.

    Methods:
        set_advanced_alignment: Set advanced alignment strategy.
        run: Execute the MUSCLE command and store the result.

    """

    def __init__(self, input_file: str, output_file: str):
        self.input_file = input_file
        self.output_file = output_file
        self.strategy: Optional[str] = None
        self.replicates: Optional[int] = None
        self.perm: Optional[str] = None
        self.perturb: Optional[int] = None
        self.results: Optional[str] = None

    def set_advanced_alignment(
        self,
        strategy: Optional[str] = None,
        replicates: Optional[int] = None,
        perm: Optional[str] = None,
        perturb: Optional[int] = None,
    ) -> None:
        """Set advanced alignment strategy."""
        if strategy:
            assert strategy in ["stratified", "diversified"], "Invalid strategy"
            assert not perm and not perturb, "Cannot set -perm or -perturb with -stratified or -diversified"
            self.strategy = strategy

        self.replicates = replicates
        self.perm = perm
        self.perturb = perturb

    def _build_command(self) -> list[str]:
        """Construct the MUSCLE command."""
        cmd = ["muscle", "-align", f"{self.input_file}", "-output", f"{self.output_file}"]

        if self.strategy:
            cmd.append(f"-{self.strategy}")
        if self.replicates:
            cmd.extend(["-replicates", f"{self.replicates}"])
        if self.perm:
            cmd.extend(["-perm", f"{self.perm}"])
        if self.perturb is not None:
            cmd.extend(["-perturb", f"{self.perturb}"])

        return cmd

    def run(self) -> None:
        """Execute the MUSCLE command and store the result."""
        cmd = self._build_command()
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            raise MuscleFailed(f"MUSCLE failed: {stderr.decode().strip()}")

        self.results = stdout.decode().strip()
