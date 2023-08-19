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
from typing import Dict, Iterator, List

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from numpy.typing import NDArray


class MSAParser:
    """A parser for multiple sequence alignment (MSA) files.

    Attributes:
        seq_ids (List[str]): The sequence IDs.
        sequences (Dict[str, Seq]): The sequences.

    """

    def __init__(self, filepath: str) -> None:
        self._sequences: Dict[str, Seq] = {}
        self._seq_ids: List[str] = []
        self._parse_file(filepath)

    def _parse_file(self, filepath: str) -> None:
        for record in SeqIO.parse(filepath, "fasta"):
            assert isinstance(record, SeqRecord)
            self._seq_ids.append(record.id)
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

    def get_seq_ids(self) -> List[str]:
        """Get the sequence IDs.

        Returns:
            List[str]: The sequence IDs.

        """
        return self._seq_ids

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
    def seq_ids(self) -> List[str]:
        """Get the sequence IDs.

        Returns:
            List[str]: The sequence IDs.

        """
        return self._seq_ids

    @property
    def sequences(self) -> Dict[str, Seq]:
        """Get the sequences.

        Returns:
            Dict[str, Seq]: The sequences.

        """
        return self._sequences

    def __getitem__(self, seq_id: str) -> Seq:
        return self._sequences[seq_id]

    def __len__(self) -> int:
        return len(self._sequences)

    def __iter__(self) -> Iterator[str]:
        return iter(self._sequences)
