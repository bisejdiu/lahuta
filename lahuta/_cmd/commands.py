from typing import Any, Mapping, get_type_hints

from typing_extensions import Required

from .base import CommandWithSubCommand, TestBaseCommand
from .options import (
    CLIOptions,
    ConvertAlisOptions,
    ConvertAlisOptionsDefaults,
    CreateDBOptions,
    CreateDBOptionsDefaults,
    SearchOptions,
    SearchOptionsDefaults,
    StructureAlignOptionsDefaults,
    StructureAlignptions,
)


def get_required_keys(typed_dict: type[CLIOptions]) -> list[str]:
    required_keys: list[str] = []  # need to use list instead of set to perserve order
    for key, value in get_type_hints(typed_dict).items():
        # Check if the type hint is wrapped with Required
        if hasattr(value, "__origin__") and value.__origin__ is Required:
            required_keys.append(key)
    return required_keys


class FoldSeekCommand(CommandWithSubCommand):
    def __init__(
        self, command_name: str, options: CLIOptions, options_defaults: dict[str, str], options_def: type
    ) -> None:
        required_keys = get_required_keys(options_def)

        # Validate required arguments
        if not set(required_keys).issubset(options):
            raise ValueError(f"Required arguments {required_keys} are missing.")

        args: list[str] = [options.pop(key) for key in required_keys]  # type: ignore
        opts: Mapping[str, Any] = {**options_defaults, **options}

        super().__init__("foldseek", command_name, args, opts)


class CreateDBCommand(FoldSeekCommand):
    def __init__(self, *, options: CreateDBOptions) -> None:
        super().__init__("createdb", options, CreateDBOptionsDefaults, CreateDBOptions)


class SearchCommand(FoldSeekCommand):
    def __init__(self, *, options: SearchOptions) -> None:
        super().__init__("search", options, SearchOptionsDefaults, SearchOptions)


class ConvertAlisCommand(FoldSeekCommand):
    def __init__(self, *, options: ConvertAlisOptions) -> None:
        super().__init__("convertalis", options, ConvertAlisOptionsDefaults, ConvertAlisOptions)


class StructureAlignCommand(FoldSeekCommand):
    def __init__(self, *, options: ConvertAlisOptions) -> None:
        super().__init__("structurealign", options, StructureAlignOptionsDefaults, StructureAlignptions)


class TestCreateDBCommand(TestBaseCommand):
    NAME = "foldseek"

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__("foldseek", "createdb", ["1gzm.pdb", "db/query"], options)


class TestSearchCommand(TestBaseCommand):
    NAME = "foldseek"

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__("foldseek", "search", ["db/query", "db/target", "db/result", "db/search_tmp"], options)


class TestConvertAlisCommand(TestBaseCommand):
    NAME = "foldseek"

    def __init__(self, options: dict[str, str]) -> None:
        super().__init__("foldseek", "convertalis", ["db/query", "db/target", "db/result", "x_aln_x_.8"], options)
