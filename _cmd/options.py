from typing import Literal
from typing_extensions import TypedDict, Required

FilePath = str | list[str]

class CLIOptions(TypedDict, total=False):
    ...

RequiredOpts = Literal["input_files"] | Literal["db_out_path"]

class CreateDBOptions(CLIOptions, total=False):
    input_files: Required[FilePath]
    db_out_path: Required[str]
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
    # 'input_files': "",
    # 'db_out_path': "",
    'chain_name_mode': "0",
    'write_mapping': "0",
    'mask_bfactor_threshold': "0.0",
    'threads': "4",
    'coord_store_mode': "2",
    'write_lookup': "1",
    'tar_include': ".*",
    'tar_exclude': "^$",
    'file_include': ".*",
    'file_exclude': "^$",
    'v': "3",
}