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
from typing import Any, Iterable, Iterator, Literal, Optional, Type, TypeVar

import numpy as np
from Bio.Seq import Seq
from numpy.typing import NDArray

from lahuta.msa.io import read_fasta
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
        self,
        *,
        filepath: Optional[str] = None,
        sequences: Optional[dict[str, Seq] | dict[str, str]] = None,
        _ref_alig_keys: Optional[set[str]] = None,
    ) -> None:
        self._labels: Optional[NDArray[np.str_]] = None
        self._sequences: dict[str, Seq] = {}
        self._ref_alig: Optional[MSAParser] = None
        match (filepath, sequences):
            case (None, _):
                assert sequences is not None
                self._sequences = {k: Seq(v) for k, v in sequences.items()}
            case (_, None):
                assert filepath is not None
                self._sequences = read_fasta(filepath)
            case (_, _):
                raise ValueError("Either filepath or sequences must be provided")

        self._is_same_length: bool = self.all_same_length(list(self._sequences.values()))

        if _ref_alig_keys is not None:
            self._ref_alig = type(self)(sequences={k: self._sequences[k] for k in _ref_alig_keys})

    # TODO (bisejdiu): rename method
    @staticmethod
    def to_indices_array(seq: str | Seq) -> NDArray[np.int32]:
        """Convert a sequence to an array of indices.

        Args:
            seq (Seq): The sequence.

        Returns:
            NDArray[np.int32]: The array of indices.

        """
        if isinstance(seq, str):
            seq = Seq(seq)
        arr = MSAParser._to_framebuffer(seq)
        return np.nonzero(arr != 45)[0]  # 45 is the ASCII code for '-'

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
        return type(self)(filepath=aligner.output_file, _ref_alig_keys=aligner.ref_alig_keys)

    @staticmethod
    def map_labels(labels: Iterable[str], sequences: list[str | Seq], fill: str = "-") -> NDArray[np.str_]:
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

    def assign_labels(self, ref_sequences: dict[str, Seq], labels: Iterable[str], fill: str = "-") -> None:
        """Assign labels to the aligned sequences.

        Args:
            ref_sequences (dict[str, Seq]): The reference sequences used for alignment.
            labels (Iterable[str]): The labels.
            fill (str, optional): The fill character. Defaults to "-".

        """
        seq_block = np.array([list(seq) for seq in ref_sequences.values()])
        ref_seq_block = np.array([list(seq) for key, seq in self.sequences.items() if key in ref_sequences])
        comp_labels = np.full(ref_seq_block.shape[1], fill, dtype="U25")

        non_fill_indices, ref_non_fill_indices = (
            np.any(seq_block != fill, axis=0),
            np.any(ref_seq_block != fill, axis=0),
        )

        non_fill_labels = np.array(labels)[non_fill_indices]
        comp_labels[ref_non_fill_indices] = non_fill_labels

        self._labels = comp_labels

    def get_labels(
        self,
        *,
        source_seq_id: Optional[str] = None,
        source_seq: Optional[str | Seq] = None,
        target_seq_ids: Optional[list[str]] = None,
        target_seqs: Optional[list[str | Seq]] = None,
        fill: str = "-",
    ) -> NDArray[np.str_]:
        """Get the labels.

        Given a sequence or sequence ID, get the labels for the residues that are not dashes in the sequence.
        The labels are taken from all existing sequences. To get the labels only from a subset of sequences,
        use the `target_seq_ids` or `target_seqs` arguments. This will return the labels for the residues that
        are not dashes in the source sequence and in any of the target sequences.

        Args:
            source_seq_id (Optional[str], optional): The source sequence ID. Defaults to None.
            source_seq (Optional[str | Seq], optional): The source sequence. Defaults to None.
            target_seq_ids (Optional[list[str]], optional): The target sequence IDs. Defaults to None.
            target_seqs (Optional[list[str | Seq]], optional): The target sequences. Defaults to None.
            fill (str, optional): The fill character. Defaults to "-".

        Returns:
            NDArray[np.str_]: The labels.

        Raises:
            ValueError: If no labels have been assigned.
            ValueError: If not all sequences are of the same length.
            ValueError: If both `source_seq_id` and `source_seq` are provided.
            ValueError: If both `target_seq_ids` and `target_seqs` are provided.
            ValueError: If the source sequence or ID is invalid.

        """
        if self._labels is None:
            raise ValueError("Cannot get labels because no labels have been assigned!")

        if not self._is_same_length:
            raise ValueError("Cannot get labels because not all sequences are of the same length!")

        self._validate_mutually_exclusive(
            source_seq_id, source_seq, "Provide either source_seq_id or source_seq, not both."
        )
        self._validate_mutually_exclusive(
            target_seq_ids, target_seqs, "Provide either target_seq_ids or target_seqs, not both."
        )

        seq = Seq(source_seq) if source_seq is not None else self._sequences.get(source_seq_id)  # type: ignore
        if seq is None:
            raise ValueError("Invalid source sequence or ID.")

        comp_seq_block = np.array(self._get_sequence_block(target_seq_ids, target_seqs))

        seq_non_fill = np.array(seq) != fill
        seq_block_non_fill = np.any(comp_seq_block != fill, axis=0)
        retain_labels_condition = seq_non_fill & seq_block_non_fill
        result_labels = np.where(retain_labels_condition, self._labels, fill)

        return result_labels

    def save(self, file_name: str) -> None:
        """Write the sequences to a FASTA file.

        Args:
            file_name (str): The file name.

        """
        with Path(file_name).open("w") as f:
            for key, value in self.sequences.items():
                f.write(f">{key.split('.')[0]}\n{value}\n")

    @staticmethod
    def all_same_length(strings: list[str | Seq]) -> bool:
        """Check if all strings are of the same length.

        Args:
            strings (list[str]): The strings.

        Returns:
            bool: True if all strings are of the same length, False otherwise.

        """
        first_length = len(strings[0])
        return all(len(s) == first_length for s in strings)

    @property
    def labels(self) -> NDArray[np.str_]:
        """Get the labels.

        Returns:
            NDArray[np.str_]: The labels.

        """
        if self._labels is None:
            raise ValueError("No labels have been assigned!")
        return self._labels

    def _validate_mutually_exclusive(self, param1: Any, param2: Any, error_message: str) -> None:  # noqa: ANN401
        if param1 is not None and param2 is not None:
            raise ValueError(error_message)

    def _get_sequence_block(self, ids: Optional[list[str]], seqs: Optional[list[str | Seq]]) -> list[list[str]]:
        if ids is not None:
            return [list(self._sequences[_id_]) for _id_ in ids]

        if seqs is not None:
            return [list(str(seq)) for seq in seqs]

        return [list(str(seq)) for seq in self._sequences.values()]

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
