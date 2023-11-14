from typing import cast, get_type_hints
from typing_extensions import Required, LiteralString

from cmd_base import Command
from options import CreateDBOptions, CreateDBOptionsDefaults, CLIOptions


def get_required_keys(typed_dict: type[CLIOptions]) -> set[LiteralString]:
    required_keys = set()
    for key, value in get_type_hints(typed_dict).items():
        # Check if the type hint is wrapped with Required
        if hasattr(value, '__origin__') and value.__origin__ is Required:
            required_keys.add(key)
    return required_keys

class CreateDBCommand(Command):
    def __init__(self, *, options: CreateDBOptions) -> None:

        required_keys = get_required_keys(CreateDBOptions)

        # Validate required arguments
        if not required_keys.issubset(options):
            raise ValueError(f"Required arguments {required_keys} are missing.")

        args, self.options = [], {**CreateDBOptionsDefaults.copy()}
        for key in CreateDBOptionsDefaults:
            if key in required_keys:
                args.append(options.get(key))
                self.options.pop(key)
            else:
                self.options[key] = cast(str, options.get(key)) or CreateDBOptionsDefaults[key]

        super().__init__("createdb", args, self.options)

createdb = CreateDBCommand(options={
    "input_files": "1gzm.pdb",
    "db_out_path": "db/query",
    # "chain_name_mode": 1,
    })

print(createdb.arguments)
