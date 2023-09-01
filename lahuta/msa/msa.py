"""MSA parser module.

This module provides a parser for multiple sequence alignment (MSA) files.
The parser is designed to work with FASTA files.

Classes:
    MSAParser: A parser for multiple sequence alignment (MSA) files.

Example:
    parser = MSAParser('path/to/msa.fasta')
    seq_ids = parser.get_seq_ids()
    seq_id = seq_ids[0]

"""
from typing import Iterator, Literal, Optional, Type, TypeVar

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from numpy.typing import NDArray

from lahuta.msa.mafft import Mafft
from lahuta.msa.muscle import Muscle

__all__ = ["MSAParser"]

T = TypeVar("T", bound="MSAParser")


MSAFactory: dict[str, Type[Mafft | Muscle]] = {
    "mafft": Mafft,
    "muscle": Muscle,
}


class InvalidBackendError(Exception):
    pass


def get_backend(backend: Literal["mafft", "muscle"] = "mafft") -> Type[Mafft | Muscle]:
    Align = MSAFactory.get(backend)
    if Align is None:
        raise InvalidBackendError(f"Invalid backend: {backend}! Supported backends: mafft, muscle")
    return Align


class MSAParser:
    """A parser for multiple sequence alignment (MSA) files.

    Attributes:
        seq_ids (list[str]): The sequence IDs.
        sequences (dict[str, Seq]): The sequences.

    """

    def __init__(
        self, filepath: Optional[str] = None, sequences: Optional[dict[str, Seq] | dict[str, str]] = None
    ) -> None:
        self._sequences: dict[str, Seq] = {}
        match (filepath, sequences):
            case (None, _):
                assert sequences is not None
                self._sequences = {k: Seq(v) for k, v in sequences.items()}
            case (_, None):
                assert filepath is not None
                self._parse_file(filepath)
            case (_, _):
                raise ValueError("Either filepath or sequences must be provided")

    def _parse_file(self, filepath: str) -> None:
        for record in SeqIO.parse(filepath, "fasta"):
            self._sequences[record.id] = record.seq

    def sequence_indices(self, seq_id: str) -> NDArray[np.int32]:
        """Get the indices of the sequence.

        Args:
            seq_id (str): The ID of the sequence.

        Returns:
            NDArray[np.int32]: The indices of the sequence.

        Raises:
            ValueError: If the sequence ID is not found.

        """
        seq = self._sequences.get(seq_id, None)
        if seq is None:
            raise ValueError(f"Sequence with ID {seq_id} not found.")
        return self.to_indices_array(seq)

    @staticmethod
    def to_indices_array(seq: Seq) -> NDArray[np.int32]:
        """Convert a sequence to an array of indices.

        Args:
            seq (Seq): The sequence.

        Returns:
            NDArray[np.int32]: The array of indices.

        """
        arr = MSAParser._to_framebuffer(seq)
        return np.nonzero(arr != 45)[0]

    @staticmethod
    def _to_framebuffer(seq: Seq) -> NDArray[np.ubyte]:
        return np.frombuffer(bytes(seq), dtype=np.ubyte)

    @property
    def seq_ids(self) -> list[str]:
        """Get the sequence IDs.

        Returns:
            list[str]: The sequence IDs.

        """
        return list(self._sequences.keys())

    @property
    def sequences(self) -> dict[str, Seq]:
        """Get the sequences.

        Returns:
            dict[str, Seq]: The sequences.

        """
        return self._sequences

    def align(
        self: T,
        backend: Literal["mafft", "muscle"] = "mafft",
        n_jobs: int = 1,
        ref_alignment: str | dict[str, str] | dict[str, Seq] | None = None,
    ) -> T:
        """Align the sequences.

        Args:
            backend (Literal["mafft", "muscle"], optional): The alignment backend. Defaults to "mafft".
            n_jobs (int, optional): The number of jobs. Defaults to 1.
            ref_alignment (str | dict[str, str] | dict[str, Seq], optional): The reference alignment. Defaults to None.

        Returns:
            T: An `MSAParser` instance with the aligned sequences.

        Raises:
            ValueError: If the backend is not supported.

        """
        Align = get_backend(backend)
        aligner = Align(self._sequences, ref_alignment=ref_alignment)
        aligner.run(n_jobs=n_jobs)
        return type(self)(aligner.output_file)

    def __add__(self: T, other: T) -> T:
        return type(self)(sequences={**self.sequences, **other.sequences})

    def __sub__(self: T, other: T) -> T:
        return type(self)(sequences={k: v for k, v in self.sequences.items() if k not in other.sequences})

    def __getitem__(self, seq_id: str) -> Seq:
        return self._sequences[seq_id]

    def __len__(self) -> int:
        return len(self._sequences)

    def __iter__(self) -> Iterator[str]:
        return iter(self._sequences)
