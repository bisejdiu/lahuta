from typing import Dict, Iterable, List

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from numpy.typing import NDArray


class MSAParser:
    def __init__(self, filepath: str) -> None:
        self._sequences: Dict[str, Seq] = {}
        self._seq_ids: List[str] = []
        self._parse_file(filepath)

    def _parse_file(self, filepath: str) -> None:
        for record in SeqIO.parse(filepath, "fasta"):  # type: ignore
            assert isinstance(record, SeqRecord)
            self._seq_ids.append(record.id)
            self._sequences[record.id] = record.seq

    def sequence_indices(self, seq_id: str) -> NDArray[np.int32]:
        seq = self._sequences.get(seq_id, None)
        if seq is None:
            raise ValueError(f"Sequence with ID {seq_id} not found.")
        return self.to_indices_array(seq)

    def get_seq_ids(self) -> List[str]:
        return self._seq_ids

    @staticmethod
    def to_indices_array(seq: Seq) -> NDArray[np.int32]:
        arr = MSAParser._to_framebuffer(seq)
        return np.nonzero(arr != 45)[0]

    @staticmethod
    def _to_framebuffer(seq: Seq) -> NDArray[np.ubyte]:
        return np.frombuffer(bytes(seq), dtype=np.ubyte)

    @property
    def seq_ids(self) -> List[str]:
        return self._seq_ids

    @property
    def sequences(self) -> Dict[str, Seq]:
        return self._sequences

    def __getitem__(self, seq_id: str) -> Seq:
        return self._sequences[seq_id]

    def __len__(self) -> int:
        return len(self._sequences)

    def __iter__(self) -> Iterable[str]:
        return iter(self._sequences)
