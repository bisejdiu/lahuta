from typing import cast, get_type_hints
from typing_extensions import Required

from base import FoldSeekBaseCommand, TestBaseCommand
from options import (
    CreateDBOptions, 
    CreateDBOptionsDefaults, 
    CLIOptions, 
    SearchOptions, 
    SearchOptionsDefaults, 
    ConvertAlisOptions,  
    ConvertAlisOptionsDefaults
)


def get_required_keys(typed_dict: type[CLIOptions]) -> list[str]:
    required_keys: list[str] = [] # need to use list instead of set to perserve order
    for key, value in get_type_hints(typed_dict).items():
        # Check if the type hint is wrapped with Required
        if hasattr(value, '__origin__') and value.__origin__ is Required:
            required_keys.append(key)
    return required_keys

class FoldSeekCommand(FoldSeekBaseCommand):
    def __init__(self, command_name: str, options: CLIOptions, options_defaults: dict[str, str], options_type: type) -> None:
        required_keys = get_required_keys(options_type)

        # Validate required arguments
        if not set(required_keys).issubset(options):
            raise ValueError(f"Required arguments {required_keys} are missing.")

        args = [options.pop(key) for key in required_keys]  # type: ignore
        opts: CLIOptions = cast(CLIOptions, {**options_defaults, **options})

        super().__init__(command_name, args, opts)


class CreateDBCommand(FoldSeekCommand):
    def __init__(self, *, options: CreateDBOptions) -> None:
        super().__init__("createdb", options, CreateDBOptionsDefaults, CreateDBOptions)

class SearchCommand(FoldSeekCommand):
    def __init__(self, *, options: SearchOptions) -> None:
        super().__init__("search", options, SearchOptionsDefaults, SearchOptions)

class ConvertAlisCommand(FoldSeekCommand):
    def __init__(self, *, options: ConvertAlisOptions) -> None:
        super().__init__("convertalis", options, ConvertAlisOptionsDefaults, ConvertAlisOptions)


class TestCreateDBCommand(TestBaseCommand):
    NAME = "foldseek"
    def __init__(self, options: dict[str, str]) -> None:
        super().__init__("createdb", ["1gzm.pdb", "db/query"], options)

class TestSearchCommand(TestBaseCommand):
    NAME = "foldseek"
    def __init__(self, options: dict[str, str]) -> None:
        super().__init__("search", ["db/query", "db/target", "db/result", "db/search_tmp"], options)

class TestConvertAlisCommand(TestBaseCommand):
    NAME = "foldseek"
    def __init__(self, options: dict[str, str]) -> None:
        super().__init__("convertalis", ["db/query", "db/target", "db/result", "x_aln_x_.8"], options)
