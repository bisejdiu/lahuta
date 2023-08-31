"""MUSCLE wrapper."""
import logging
import os
import subprocess
from tempfile import NamedTemporaryFile
from typing import Dict, List, Optional

from Bio.Seq import Seq

__all__ = ["Muscle"]

logging.basicConfig(level=logging.INFO)


class MuscleFailed(Exception):
    """Exception raised when MUSCLE fails to align the sequences."""


class Muscle:
    """Wrapper for MUSCLE.

    Attributes:
        input_file (str): Path to the input file.
        result (str): Result of the alignment.
        strategy (str): Strategy to use for advanced alignment.
        options (dict[str, str] | dict[str, Seq]): Options to pass to MUSCLE.
        output_file (str): Path to the output file.

    Methods:
        set_option: Set an option to pass to MUSCLE.
        set_advanced_alignment: Set advanced alignment options.
        run: Run MUSCLE.
        cleanup: Delete the input and output files.
    """

    def __init__(self, sequence_data: str | dict[str, str] | dict[str, Seq]) -> None:
        match sequence_data:
            case str(path):
                self.input_file = path
            case dict(data):
                with NamedTemporaryFile(delete=False) as temp:
                    self.input_file = temp.name
                    temp.write(self.dict_to_fasta(data).encode())

        self.result: Optional[str] = None
        self.strategy: Optional[str] = None
        self.options: Dict[str, str] = {}
        self.output_file = self.input_file + "_aligned.fasta"

    def set_option(self, key: str, value: Optional[str]) -> None:
        """Set an option to pass to MUSCLE.

        Args:
            key (str): Option name.
            value (Optional[str]): Option value.

        """
        if value is not None:
            self.options[key] = str(value)

    def set_advanced_alignment(
        self,
        strategy: Optional[str] = None,
        replicates: Optional[int] = None,
        perm: Optional[str] = None,
        perturb: Optional[int] = None,
    ) -> None:
        """Set advanced alignment options.

        Args:
            strategy (Optional[str]): Strategy to use for advanced alignment.
            replicates (Optional[int]): Number of replicates.
            perm (Optional[str]): Permutation seed.
            perturb (Optional[int]): Perturbation seed.

        """
        if strategy:
            assert strategy in ["stratified", "diversified"], "Invalid strategy"
            assert not perm and not perturb, "Cannot set -perm or -perturb with -stratified or -diversified"
            self.strategy = strategy

        if replicates is not None:
            self.options["replicates"] = str(replicates)
        if perm is not None:
            self.options["perm"] = perm
        if perturb is not None:
            self.options["perturb"] = str(perturb)

    def _build_command(self) -> List[str]:
        cmd = ["muscle", "-align", self.input_file, "-output", self.output_file]
        if self.strategy:
            cmd.append(f"-{self.strategy}")
        for k, v in self.options.items():
            cmd.extend([f"-{k}", v])
        return cmd

    def run(self) -> None:
        """Run MUSCLE."""
        cmd = self._build_command()
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, stderr = process.communicate()

        if process.returncode != 0:
            raise MuscleFailed(f"MUSCLE failed: {stderr.decode().strip()}")

        logging.info(f"Alignment completed. Output file located at: {self.output_file}")

        with open(self.output_file, "r") as f:
            self.result = f.read()

    @staticmethod
    def dict_to_fasta(seq_dict: dict[str, str] | dict[str, Seq]) -> str:
        """Convert a dictionary of sequences to a FASTA string.

        Args:
            seq_dict (dict[str, str] | dict[str, Seq]): Dictionary of sequences.

        Returns:
            str: FASTA string.

        """
        return "".join([f">{k}\n{v}\n" for k, v in seq_dict.items()])

    def cleanup(self) -> None:
        """Delete the input and output files."""
        if os.path.exists(self.input_file):
            os.remove(self.input_file)
        if os.path.exists(self.output_file):
            os.remove(self.output_file)

    # seems to not be recommended
