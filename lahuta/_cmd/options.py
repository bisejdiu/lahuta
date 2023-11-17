from enum import Enum, auto
from typing import Literal

from typing_extensions import Required, TypedDict

FilePath = str | list[str]


class CLIOptions(TypedDict, total=False):
    ...


class CreateDBOptions(CLIOptions, total=False):
    input_files: Required[FilePath]
    db_out_path: Required[str]
    # optionals
    chain_name_mode: Literal["0", "1"]
    write_mapping: Literal["0", "1"]
    mask_bfactor_threshold: str
    threads: str
    coord_store_mode: Literal["1", "2"]
    write_lookup: Literal["0", "1"]
    tar_include: str
    tar_exclude: str
    file_include: str
    file_exclude: str
    v: Literal["0", "1", "2", "3"]


CreateDBOptionsDefaults: dict[str, str] = {
    "chain_name_mode": "0",
    "write_mapping": "0",
    "mask_bfactor_threshold": "0",
    "threads": "4",
    "coord_store_mode": "2",
    "write_lookup": "1",
    "tar_include": ".*",
    "tar_exclude": "^$",
    "file_include": ".*",
    "file_exclude": "^$",
    "v": "3",
}


class SearchOptions(CLIOptions, total=False):
    query: Required[str]
    target: Required[str]
    result: Required[str]
    search_tmp: Required[str]
    a: Literal["0", "1"]
    # optionals
    alignment_mode: Literal["0", "1", "2", "3"]
    alignment_output_mode: Literal["0", "1", "2", "3", "4", "5"]
    comp_bias_corr: str
    gap_open: str
    gap_extend: str
    s: str
    k: str
    mask: Literal["0", "1"]
    mask_prob: str
    remove_tmp_files: Literal["0", "1"]
    max_seqs: str
    exhaustive_search: Literal["0", "1"]


SearchOptionsDefaults: dict[str, str] = {
    "a": "0",
    "alignment_mode": "3",
    "alignment_output_mode": "0",
    "comp_bias_corr": "1",
    "gap_open": "aa:10,nucl:10",
    "gap_extend": "aa:1,nucl:1",
    "s": "9.5",
    "k": "0",
    "mask": "0",
    "mask_prob": "1.000",
    "remove_tmp_files": "1",
    "max_seqs": "1000",
    "exhaustive_search": "0",
}


class ConvertAlisOptions(CLIOptions, total=False):
    query: Required[str]
    target: Required[str]
    result: Required[str]
    output: Required[str]
    # optionals
    sub_mat: str
    format_mode: Literal["0", "1", "2", "3", "4", "5"]
    format_output: str
    translation_table: str
    gap_open: str
    gap_extend: str
    db_output: Literal["0", "1"]
    db_load_mode: Literal["0", "1"]
    search_type: Literal["0", "1"]
    threads: str
    compressed: Literal["0", "1"]
    v: Literal["0", "1", "2", "3"]


class FormatOutput(Enum):
    QUERY = auto()
    TARGET = auto()
    EVALUE = auto()
    GAPOPEN = auto()
    PIDENT = auto()
    FIDENT = auto()
    NIDENT = auto()
    QSTART = auto()
    QEND = auto()
    QLEN = auto()
    TSTART = auto()
    TEND = auto()
    TLEN = auto()
    ALNLEN = auto()
    RAW = auto()
    BITS = auto()
    CIGAR = auto()
    QSEQ = auto()
    TSEQ = auto()
    QHEADER = auto()
    THEADER = auto()
    QALN = auto()
    TALN = auto()
    MISMATCH = auto()
    QCOV = auto()
    TCOV = auto()
    QSET = auto()
    QSETID = auto()
    TSET = auto()
    TSETID = auto()
    TAXID = auto()
    TAXNAME = auto()
    TAXLINEAGE = auto()
    LDDT = auto()
    LDDTFULL = auto()
    QCA = auto()
    TCA = auto()
    T = auto()
    U = auto()
    QTMSCORE = auto()
    TTMSCORE = auto()
    ALNTMSCORE = auto()
    RMSD = auto()
    PROB = auto()
    QCOMPLEXTMSCORE = auto()
    TCOMPLEXTMSCORE = auto()
    ASSIGNID = auto()


class FormatOutputController:
    DEFAULTS = [
        FormatOutput.QUERY,
        FormatOutput.TARGET,
        FormatOutput.FIDENT,
        FormatOutput.ALNLEN,
        FormatOutput.MISMATCH,
        FormatOutput.GAPOPEN,
        FormatOutput.QSTART,
        FormatOutput.QEND,
        FormatOutput.TSTART,
        FormatOutput.TEND,
        FormatOutput.EVALUE,
        FormatOutput.BITS,
    ]

    def __init__(self):
        self.state = {option: False for option in FormatOutput}
        self.state.update({option: True for option in self.DEFAULTS})

    def set_on(self, options: FormatOutput) -> None:
        self.state[options] = True

    def set_off(self, options: FormatOutput) -> None:
        self.state[options] = False

    def set_all_on(self) -> None:
        for options in FormatOutput:
            self.state[options] = True

    def set_all_off(self) -> None:
        for options in FormatOutput:
            self.state[options] = False

    def get_on(self) -> str:
        # Join the names of the colors that are 'on'
        return ",".join(option.name.lower() for option, is_on in self.state.items() if is_on)

    def get_off(self) -> str:
        # Join the names of the colors that are 'off'
        return ",".join(option.name.lower() for option, is_on in self.state.items() if not is_on)

    def get_state(self) -> dict[FormatOutput, bool]:
        return self.state


ConvertAlisOptionsDefaults: dict[str, str] = {
    "sub_mat": "aa:3di.out,nucl:3di.out",
    "format_mode": "0",
    "format_output": FormatOutputController().get_on(),
    "translation_table": "1",
    "gap_open": "aa:10,nucl:10",
    "gap_extend": "aa:1,nucl:1",
    "db_output": "0",
    "db_load_mode": "0",
    "search_type": "0",
    "threads": "4",
    "compressed": "0",
    "v": "3",
}
