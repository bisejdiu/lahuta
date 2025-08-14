"""Wrapper for MAFFT."""
import logging
import subprocess
from tempfile import NamedTemporaryFile
from typing import Optional

from Bio.Seq import Seq

from lahuta.msa.io import read_keys_from_fasta

__all__ = ["Mafft"]

logging.basicConfig(level=logging.INFO)


class Mafft:
    """Wrapper for MAFFT.

    Attributes:
        input_file (str): Path to the input file.
        result (str): Result of the alignment.
        options (dict[str, str | None]): Options to pass to MAFFT.
        output_file (str): Path to the output file.

    Methods:
        set_advanced_options: Set advanced alignment options.
        run: Run MAFFT.
        read_result: Read the result from the output file.
    """

    def __init__(
        self,
        sequence_data: str | dict[str, str] | dict[str, Seq],
        ref_alignment: Optional[str | dict[str, str] | dict[str, Seq]] = None,
    ) -> None:
        self.ref_alig_keys: Optional[set[str]] = None
        match sequence_data:
            case str(path):
                self.input_file = path
            case dict(data):
                with NamedTemporaryFile(delete=False) as temp:
                    self.input_file = temp.name
                    temp.write(self.dict_to_fasta(data).encode())
            case _:
                raise ValueError("Invalid sequence data")

        match ref_alignment:
            case str(path):
                ref_alignment = path
                self.ref_alig_keys = read_keys_from_fasta(path)
            case dict(data):
                with NamedTemporaryFile(delete=False) as temp:
                    ref_alignment = temp.name
                    temp.write(self.dict_to_fasta(data).encode())
                    self.ref_alig_keys = set(data.keys())
            case _:
                ref_alignment = None

        self.ref_alignment = ref_alignment

        self._result: Optional[str] = None
        self.options: dict[str, str | None] = {}
        self._command: list[str] = []
        with NamedTemporaryFile(delete=False) as temp:
            self.output_file = temp.name

    @staticmethod
    def dict_to_fasta(data: dict[str, str] | dict[str, Seq]) -> str:
        """Convert a dictionary to a FASTA string.

        Args:
            data (dict[str, str]): Dictionary of sequences.

        Returns:
            str: FASTA string.

        """
        return "\n".join([f">{key}\n{value}" for key, value in data.items()])

    def set_advanced_options(self, *args: str, **kwargs: str | None) -> None:
        """Set advanced alignment options.

        Args:
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.

        """
        if args:
            for arg in args:
                self.options[arg] = None
        self.options.update(kwargs)

    def _build_command(self) -> list[str]:
        command = ["mafft"]

        if self.options:
            for key, value in self.options.items():
                command.append(f"--{key}")
                if value is not None:
                    command.append(str(value))

        command.append(self.input_file)

        return command

    def run(self, n_jobs: int = 1, keeplength: bool = False) -> None:
        """Run MAFFT.

        Args:
            n_jobs (int): Number of workers.
            keeplength (bool): Keep the original sequence lengths.

        """
        match (self.ref_alignment, keeplength):
            case (None, _):
                command = self._build_command()
            case (_, True):
                assert self.ref_alignment is not None
                command = ["mafft", "--add", self.input_file, "--keeplength", self.ref_alignment]
            case (_, False):
                assert self.ref_alignment is not None
                command = ["mafft", "--add", self.input_file, self.ref_alignment]
            case _:
                raise ValueError("Invalid combination of arguments")

        command[1:1] = ["--thread", str(n_jobs)]

        self._command = command
        with open(self.output_file, "w") as out_file:
            process = subprocess.Popen(command, stdout=out_file, stderr=subprocess.PIPE)
            _, logs = process.communicate()

        self.logs = logs.decode()

        logging.info(f"Alignment completed. Output file located at: {self.output_file}")

    @property
    def command(self) -> str:
        """Get the command.

        Returns:
            str: The command.

        """
        return " ".join(self._command)

    @property
    def result(self) -> str:
        """Read the result from the output file.

        Returns:
            str: Result of the alignment.

        """
        with open(self.output_file, "r") as file:
            self._result = file.read()
        return self._result
