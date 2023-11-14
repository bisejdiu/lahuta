from typing import cast, get_type_hints
from typing_extensions import Required

from cmd_base import FoldSeekCommand
from options import CreateDBOptions, CreateDBOptionsDefaults, CLIOptions


def get_required_keys(typed_dict: type[CLIOptions]) -> list[str]:
    required_keys: list[str] = [] # need to use list instead of set to perserve order
    for key, value in get_type_hints(typed_dict).items():
        # Check if the type hint is wrapped with Required
        if hasattr(value, '__origin__') and value.__origin__ is Required:
            required_keys.append(key)
    return required_keys

class CreateDBCommand(FoldSeekCommand):
    def __init__(self, *, options: CreateDBOptions) -> None:

        required_keys = get_required_keys(CreateDBOptions)

        # Validate required arguments
        if not set(required_keys).issubset(options):
            raise ValueError(f"Required arguments {required_keys} are missing.")
        
        args = [options.pop(key) for key in required_keys]  # type: ignore
        opts: CLIOptions = cast(CLIOptions, {**CreateDBOptionsDefaults, **options})

        super().__init__("createdb", args, opts)
