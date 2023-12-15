"""IO functions for MSA."""
from Bio import SeqIO
from Bio.Seq import Seq

__all__ = ["read_fasta", "read_keys_from_fasta"]

def read_fasta(path: str) -> dict[str, Seq]:
    """Read a FASTA file.

    Args:
        path (str): Path to the FASTA file.

    Returns:
        dict[str, Seq]: Dictionary of sequences.
    """
    return {record.id: record.seq for record in SeqIO.parse(path, "fasta")}

def read_keys_from_fasta(path: str) -> set[str]:
    """Read the keys from a FASTA file.

    Args:
        path (str): Path to the FASTA file.

    Returns:
        set[str]: Set of keys.
    """
    return {record.id for record in SeqIO.parse(path, "fasta")}