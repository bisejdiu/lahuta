from pathlib import Path

import pytest

import numpy as np
from Bio.Seq import Seq
from numpy.typing import NDArray

from lahuta import Luni
from lahuta.api import CachedFileProcessor, intersection, union
from lahuta.core.neighbors import LabeledNeighborPairs, NeighborPairs
from lahuta.msa import MSAParser
from lahuta.tests.base import BaseFile

single_letter_code = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}
get_single_letter_code = np.vectorize(lambda x: single_letter_code.get(x, "-"))

class R1(BaseFile):
    FILE_NAME = "1F88"
    def __init__(self, pdb: bool = True) -> None:
        super().__init__(pdb=pdb)

class R2(BaseFile):
    FILE_NAME = "1HZX"
    def __init__(self, pdb: bool = True) -> None:
        super().__init__(pdb=pdb)

TEST_PARAMS = [(R1, "A"), (R2, "A")]
MSA_PDB_FILES = [
    ("1f88", "A", 6080, 841),
    ("1hzx", "A", 6988, 893),
    ("1l9h", "A", 7299, 944),
    ("1gzm", "A", 7127, 993),
    ("1u19", "A", 7393, 1085),
    ("2hpy", "B", 6815, 998),
    ("2g87", "A", 7377, 1079),
    ("2i35", "A", 5569, 734),
    ("2i36", "A", 4340, 602),
    ("2i37", "A", 3935, 566),
]
MSA_PDB_DICT = {key: (v1, v2, v3) for key, v1, v2, v3 in MSA_PDB_FILES}

def read_luni(pdb_file: str) -> Luni:
    """Read a PDB file and return a Luni object."""
    return Luni(pdb_file)

def load_labels(non_dash_columns: NDArray[np.bool_]) -> list[str]:
    """Load labels from a file."""
    file_loc = Path(__file__).parent / "data/labels.txt"
    with open(file_loc.as_posix(), "r") as f:
        labels = f.read().splitlines()

    labels = np.array(labels)[np.where(non_dash_columns == True)[0]].tolist()
    return labels

def load_msa(file_path: str) -> MSAParser:
    """Load an MSA from a file."""
    return MSAParser(file_path)

def load_rhod_msa() -> MSAParser:
    """Load an MSA from a file."""
    file_loc = Path(__file__).parent / "data/rhodopsins.fasta"
    return MSAParser(file_loc.as_posix())

def process_neighbors(file_path: str) -> NeighborPairs:
    """Extracts the neighbors from a PDB file."""
    basename = Path(file_path).name
    chain_id = MSA_PDB_DICT[basename[:-4]][0]

    luni = Luni(file_path).filter(f"chainID {chain_id}")

    ns = luni.neighbors(radius=5.0, res_dif=4)

    return ns

def process_sequence(file_path: str) -> str:
    """Extracts the sequence from a PDB file."""
    base_name = Path(file_path).name
    chain_id = MSA_PDB_DICT[base_name[:-4]][0]
    luni = Luni(file_path).filter(f"chainID {chain_id}")
    return luni.sequence(use_synonyms=True)

@pytest.mark.parametrize("pdb_class, chain", TEST_PARAMS)
def test_read_pdb(pdb_class: type[BaseFile], chain: str) -> None:
    luni = read_luni(pdb_class().file_loc)
    assert luni.n_atoms is not None

    filtered_luni = luni.filter(f"chainID '{chain}'")
    assert filtered_luni.n_atoms is not None

def test_process_sequence() -> None:
    R1_seq = (
        "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLR"
        "TPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVV"
        "VCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYM"
        "FVVHFIIPLIVIFFCYGQLVFTVKEAAASATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQG"
        "SDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPSTTVSKTETSQVAPA"
    )
    R2_seq = (
        "XMNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKK"
        "LRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIER"
        "YVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESF"
        "VIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIF"
        "THQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDSTTVSKTETSQVAPA"
    )
    assert process_sequence(R1().file_loc) == R1_seq
    assert process_sequence(R2().file_loc) == R2_seq

    sequence_processor = CachedFileProcessor(
        file_list=[R1().file_loc, R2().file_loc],
        worker=process_sequence,
    )
    sequence_processor.process()
    assert sequence_processor.results == {
        Path(R1().file_loc).name: R1_seq,
        Path(R2().file_loc).name: R2_seq,
    }


def correct_seqs(ref_parser: MSAParser) -> tuple[MSAParser, NDArray[np.bool_]]:
    """Corrects the sequences in the reference MSA."""
    sequences = list(ref_parser.sequences.values())
    sequences_block = np.array([list(seq) for seq in sequences])

    # Find columns with at least one non-dash element, and extract their indices.
    non_dash_columns = np.any(sequences_block != "-", axis=0)

    short_ref_seqs = {}
    for seq in ref_parser.sequences:
        non_dash_seq = np.array(ref_parser.sequences[seq])[non_dash_columns]
        short_ref_seqs[seq] = "".join(non_dash_seq)

    return MSAParser(sequences=short_ref_seqs), non_dash_columns

def get_alig_seq_labels(ref_seq: str | Seq, target_seq: str | Seq, mapped_labels: NDArray[np.str_]) -> NDArray[np.str_]:

    bov_pdb = MSAParser.to_indices_array(ref_seq)
    bov_seq = MSAParser.to_indices_array(target_seq)

    intersect_indices = np.intersect1d(bov_pdb, bov_seq)

    alig_seq_labels = np.full(len(ref_seq), "-", dtype="U25")
    alig_seq_labels[intersect_indices] = mapped_labels[intersect_indices]

    return alig_seq_labels


def test_get_aligned_seqs() -> None:
    objects_store = []
    for file_name, _, _, _ in MSA_PDB_FILES:
        NewClass = type(file_name, (BaseFile,), {"FILE_NAME": file_name})
        # Instantiating the object
        obj = NewClass(pdb=True)
        objects_store.append(obj)

    # load all structures from MSA_PDB_FILES using type()
    sequence_processor = CachedFileProcessor(
        file_list=[obj.file_loc for obj in objects_store],
        worker=process_sequence,
    )
    sequence_processor.process()
    parser = MSAParser(sequences=sequence_processor.results)
    assert sorted(parser.seq_ids) == sorted([f"{x[0]}".lower() + ".pdb" for x in MSA_PDB_FILES])
    
    unprocessed_ref_parser = load_rhod_msa()
    assert len(unprocessed_ref_parser.seq_ids) == 4
    assert set(unprocessed_ref_parser.sequences.keys()) == {
        "Rhodopsin_H._adansoni",
        "Rhodopsin_Bovine",
        "Rhodopsin_Human",
        "Rhodopsin_J._flying_s.."
    }
    assert len(next(iter(unprocessed_ref_parser.sequences.values()))) == 455

    ref_parser, non_dash_columns = correct_seqs(unprocessed_ref_parser)
    ndash_bytes = (
        b"\x00\x00\x00\x0f\xff\xdf\xff\xfc\x0f\x03\xff\xff\xff\xff\x00\xd8\x00\x00"
        b"\x7f\xff\xff\xff\xf8\x01\xfe\x03\xff\xdf\xff\x81\xc0\x1f\xf5\xff\xff\xff"
        b"\xff\xfc\x01\xff\xff\xff\xfb\xff\xf8\x00\x00\x00\x1f\xff\xbf\xde\x00\x07"
        b"\xff\xc0\x00"
    )
    unpacked = np.unpackbits(np.frombuffer(ndash_bytes, dtype=np.uint8))
    unpacked = unpacked[:455]

    assert len(ref_parser.seq_ids) == 4
    assert set(ref_parser.sequences.keys()) == {
        "Rhodopsin_H._adansoni",
        "Rhodopsin_Bovine",
        "Rhodopsin_Human",
        "Rhodopsin_J._flying_s.."
    }
    assert len(next(iter(ref_parser.sequences.values()))) == 281
    assert np.all(unpacked == non_dash_columns)

    # Align the sequences
    aligned_seqs = parser.align(n_jobs=4, ref_alignment=ref_parser.sequences)
    assert set(aligned_seqs.seq_ids) == set(ref_parser.seq_ids) | set(parser.seq_ids)

    # Get labels
    labels = load_labels(non_dash_columns)
    assert len(labels) == 281

    aligned_ref_seqs = (aligned_seqs - parser)
    mapped_labels = MSAParser.map_labels(labels, list(aligned_ref_seqs.sequences.values()))

    processor = CachedFileProcessor(
        file_list=[obj.file_loc for obj in objects_store],
        worker=process_neighbors,
    )
    processor.process(n_jobs=4)

    for result, obj in processor.results.items():
        assert isinstance(obj, NeighborPairs)
        
        key = result[:-4]
        assert key in MSA_PDB_DICT
        assert obj.pairs.shape[0] == MSA_PDB_DICT[key][1]

    mapped_results = {}
    for file_name in processor.results:
        basename = Path(file_name).name
        mapped_results[basename] = processor.results[basename].map(aligned_seqs.sequences[basename], fields={"names": False, "resnames": False, "chainids": False})

    for result, obj in mapped_results.items():
        assert isinstance(obj, LabeledNeighborPairs)
        assert obj.pairs.shape[0] == MSA_PDB_DICT[result[:-4]][2]

    assert union(*mapped_results.values()).pairs.shape[0] == 2217
    assert intersection(*mapped_results.values()).pairs.shape[0] == 426 


    ref_seq = aligned_seqs.sequences["1f88.pdb"]
    target_seq = aligned_seqs.sequences["Rhodopsin_Bovine"]

    alig_seq_labels = get_alig_seq_labels(ref_seq, target_seq, mapped_labels)
    assert alig_seq_labels.shape == (366,)

    luni = Luni(R1().file_loc).filter("chainID A") # 1GZM:A
    ns = luni.neighbors(radius=5, res_dif=4)
    mapped_ns = ns.map(ref_seq, cusotm_fields={
        "gns": {
            "values": alig_seq_labels.astype(str), # type: ignore
            "fill": "-"
        },
    })

    assert mapped_ns.pairs.shape == (6080 , 2)
    is_correct_gns_mapping(mapped_ns, ref_seq, alig_seq_labels, aligned_seqs)

    for obj in objects_store:
        key = obj.FILE_NAME[:-4]
        if key in {"2hpy", "2i37"}: continue
        luni = Luni(obj.file_loc).filter(f"chainID {MSA_PDB_DICT[key][0]}")
        ns = luni.neighbors(radius=5, res_dif=4)

        ref_seq = aligned_seqs.sequences[f"{obj.FILE_NAME.lower()}"]
        alig_seq_labels = get_alig_seq_labels(ref_seq, target_seq, mapped_labels)
        mapped_ns = ns.map(ref_seq, cusotm_fields={
            "gns": {
                "values": alig_seq_labels.astype(str), # type: ignore
                "fill": "-", 
                "pdb_id": key,
            },
        })

        assert mapped_ns.pairs.shape == (MSA_PDB_DICT[key][1] , 2)
        is_correct_gns_mapping(mapped_ns, ref_seq, alig_seq_labels, aligned_seqs)


def is_correct_gns_mapping(mapped_ns: LabeledNeighborPairs, ref_seq: str | Seq, alig_seq_labels: NDArray[np.str_], aligned_seqs: MSAParser) -> None:
    arr = LabeledNeighborPairs.remove_duplicates(mapped_ns.pairs[["resids", "resnames", "gns"]])
    resids_as_int = arr["resids"].astype(int)
    seq_extended = ref_seq + "-" * (resids_as_int.max() - len(ref_seq) + 1)
    seq_extended = np.array(seq_extended)
    target_arr = seq_extended[resids_as_int]

    resnames = arr["resnames"]
    ref_array = get_single_letter_code(resnames)
    assert np.all(ref_array == target_arr)

    labels = np.full(seq_extended.size, "-", dtype="<U25")
    labels[aligned_seqs.to_indices_array(ref_seq)] = alig_seq_labels[aligned_seqs.to_indices_array(ref_seq)]
    target_arr = labels[resids_as_int]
    assert np.all(target_arr == arr["gns"])
