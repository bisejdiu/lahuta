import asyncio
import os
from pathlib import Path

from .base import CommandWithoutSubCommand, SimpleCommand
from .commands import ConvertAlisCommand, CreateDBCommand, SearchCommand
from .runner import WorkflowRunner


class EasySearchWorkflow:
    def __init__(self, query: str, target: str) -> None:
        self.query = query
        self.target = target
        self.runner = WorkflowRunner()

    def _main_command_loop(self) -> None:
        create_query_db = CreateDBCommand(
            options={
                "input_files": self.query,
                "db_out_path": "db/query",
            }
        )
        create_target_db = CreateDBCommand(
            options={
                "input_files": self.target,
                "db_out_path": "db/target",
            }
        )
        search = SearchCommand(
            options={
                "query": "db/query",
                "target": "db/target",
                "result": "db/result",
                "search_tmp": "db/search_tmp",
                "a": "1",
                "k": "6",
                "mask_prob": "0.99995",
            }
        )
        convert_alis_fm0 = ConvertAlisCommand(
            options={
                "query": "db/query",
                "target": "db/target",
                "result": "db/result",
                "output": "x_aln_x_.8",
                "format_mode": "0",
            }
        )
        convert_alis_fm5 = ConvertAlisCommand(
            options={
                "query": "db/query",
                "target": "db/target",
                "result": "db/result",
                "output": "x_aln_x_.8",
                "format_mode": "5",
            }
        )

        if not os.path.exists("db"):  # noqa: PTH110
            self.runner.add_command(CommandWithoutSubCommand("mkdir", ["db"]))
        self.runner.add_command(create_query_db)
        self.runner.add_command(create_target_db, dependencies=[create_query_db])
        self.runner.add_command(search, dependencies=[create_target_db])
        self.runner.add_command(convert_alis_fm0, dependencies=[search])
        self.runner.add_command(convert_alis_fm5, dependencies=[search])

    def seek(self) -> None:
        self._main_command_loop()
        self.state = self.runner.run()

    @property
    def output(self) -> str:
        if not self.runner.is_successful():
            return f"Workflow failed with error: \n{self.runner.error}"

        try:
            results = ""
            for output in self.state.result().values():
                results += output + "\n"
        except asyncio.InvalidStateError:
            return "Workflow is still running. Please wait..."

        return results


class SimpleBashCommandsWorkflow:
    def __init__(self, path: str | Path):
        self.path = Path(path).resolve()
        self.old_path = Path.cwd()
        self.runner = WorkflowRunner()

    def _main_command_loop(self) -> None:
        self.runner.add_command(SimpleCommand("ls"))
        self.runner.add_command(SimpleCommand("pwd"))
        self.runner.add_command(SimpleCommand("ls", {"l": ""}))

    def seek(self) -> None:
        self._main_command_loop()
        self.state = self.runner.run()

    @property
    def output(self) -> str:
        if not self.runner.is_successful():
            return f"Workflow failed with error: \n{self.runner.error}"

        try:
            results = ""
            for output in self.state.result().values():
                results += output + "\n"

        except asyncio.InvalidStateError:
            return "Workflow is still running. Please wait..."

        return results


class ChangeDirContextManager:
    def __init__(self, path: str | Path) -> None:
        self.path = Path(path).resolve()
        self.old_path = Path.cwd()

    def __enter__(self) -> None:
        os.chdir(self.path)

    def __exit__(self, exc_type: object, exc_value: object, traceback: object) -> None:
        os.chdir(self.old_path)


if __name__ == "__main__":
    # from commands import TestCreateDBCommand, TestSearchCommand, TestConvertAlisCommand
    # runner = WorkflowRunner()
    # runner.add_command(TestCreateDBCommand(options={
    #         "input_files": "query",
    #         "db_out_path": "db/query",
    #     }))
    # runner.add_command(TestSearchCommand(options={
    #         "query": "db/query",
    #         "target": "db/target",
    #         "result": "db/result",
    #         "search_tmp": "db/search_tmp"
    #     }))
    # runner.add_command(TestConvertAlisCommand(options={
    #         "query": "db/query",
    #         "target": "db/target",
    #         "result": "db/result",
    #         "output": "x_aln_x_.8",
    #         "format_mode": "0",
    #     }))

    # runner.run()

    # assert runner.is_successful() is True
    # assert list(runner.get_state().values()) == ["1", "2", "3"]

    # w = EasySearchWorkflow(query='1gzm.pdb', target='examples/')
    # w.seek()
    # print(w.output)

    # w = SimpleBashCommandsWorkflow('/mnt/f/foldseek/lahuta_tests/test_dir_change')
    # w.seek()
    # print(w.output)

    # with ChangeDirContextManager('/mnt/f/foldseek/lahuta_tests/test_dir_change'):
    #     w = SimpleBashCommandsWorkflow('/mnt/f/foldseek/lahuta_tests/test_dir_change')
    #     w.seek()
    #     print(w.output)

    with ChangeDirContextManager("/mnt/f/foldseek/lahuta_tests/test_dir_change"):
        w = EasySearchWorkflow(query="../1gzm.pdb", target="../examples/")
        w.seek()
        print(w.output)

    # runner.add_command(FoldSeekCommand("easy-search", [""], {"h": ""}))
