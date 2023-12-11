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
from pathlib import Path
from typing import Iterable, Iterator, Literal, Optional, Type, TypeVar

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

    @sequences.setter
    def sequences(self, sequences: dict[str, Seq]) -> None:
        """Set the sequences.

        Args:
            sequences (dict[str, Seq]): The sequences.

        """
        self._sequences = sequences

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

    @staticmethod
    def map_labels(labels: Iterable[str], sequences: list[str | Seq], fill: str="-") -> NDArray[np.str_]:
        """Map labels to the aligned reference sequences.

        Identifies the unique indices that have been mapped in the sequences, creates an array of the same length 
        as the reference sequence filled with the `fill` value, and then places the labels at the appropriate 
        mapped indices.

        Args:
            labels (Iterable[str]): The labels to be mapped. The size of labels should match the size of unaligned 
                reference sequences.
            sequences (list[str]): A list of aligned reference sequences. All sequences must be of 
                the same length.
            fill (str, optional): The character used to fill unmapped positions in the output array. Defaults to '-'.

        Returns:
            NDArray[str]: An array of the same length as the reference sequences, where each position is either filled 
                with a label (if that index is mapped) or the specified fill character.

        Raises:
            AssertionError: If the aligned reference sequences are not of the same length.

        Example:
            labels = ['1', '2', '3', '4']
            sequences = [Seq('--ABC'), Seq('A-BC-')]
            mapped_labels = map_labels(labels, sequences)
            # mapped_labels would be ['1', '-', '2', '3', '4']
        """
        # all sequences should have the same length
        assert len({len(seq) for seq in sequences}) == 1
        
        # Find all the unique indices that have been mapped.
        mapped_indices = np.concatenate([MSAParser.to_indices_array(seq) for seq in sequences])
        mapped_indices = np.unique(mapped_indices)

        # Map the labels to the indices.
        mapped_labels = np.full(len(sequences[0]), fill, dtype=object)
        mapped_labels[mapped_indices] = labels
        return mapped_labels.astype(str)

    @staticmethod
    def map_labels_alt(labels: Iterable[str], sequences: list[str | Seq], fill: str="-") -> NDArray[np.str_]:
        """Map labels to positions in sequences where at least one sequence differs from the fill character.

        Identifies columns in the sequence list where at least one element is not the fill character,
        and maps the provided labels to these columns. It creates an array, the length of a sequence, initially filled
        with the fill character, and replaces the fill character at these specific columns with the labels.

        Args:
            labels (Iterable[str]): The labels to be mapped. The number of labels should match the number of 
                non-fill positions (non-dash columns) in the sequences.
            sequences (list[str]): A list of aligned sequences where the fill character represents
                gaps or unmapped positions.
            fill (str, optional): The character used to identify unmapped or gap positions in the sequences. 
                Defaults to '-'.

        Returns:
            NDArray[np.str_]: An array of the same length as the sequences, where positions corresponding to non-fill 
                columns in the sequences are replaced with the provided labels, and other positions are filled 
                with the specified fill character.

        Example:
            labels = ['1', '2', '3', '4']
            aligned_ref_sequences = [Seq('--ABC'), Seq('A-BC-')]
            mapped_labels = map_labels(labels, aligned_ref_sequences)
            # mapped_labels would be ['1', '-', '2', '3', '4']

        Note:
            The function assumes that all sequences in the list are of the same length. It does not perform 
            any alignment but relies on the sequences already being aligned, with the fill character indicating 
            gaps or unmapped positions.
        """
        sequences_block = np.array([list(seq) for seq in sequences])

        # Find columns with at least one non-dash element, and extract their indices.
        non_dash_columns = np.any(sequences_block != fill, axis=0)

        # Map the labels to the indices.
        mapped_labels = np.full(sequences_block.shape[1], fill, dtype=object)
        mapped_labels[non_dash_columns] = labels
        return mapped_labels

    def save(self, file_name: str) -> None:
        """Write the sequences to a FASTA file.

        Args:
            file_name (str): The file name.

        """
        with Path(file_name).open("w") as f:
            for key, value in self.sequences.items():
                f.write(f">{key.split('.')[0]}\n{value}\n")

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
