from asyncer import WorkflowRunner
from base import BaseCommand, FoldSeekBaseCommand
from createdb import CreateDBCommand, SearchCommand, ConvertAlisCommand
from createdb import TestCreateDBCommand, TestSearchCommand, TestConvertAlisCommand


class EasySearchWorkflow:
    def __init__(self, query, target):
        self.query = query
        self.target = target
        self.runner = WorkflowRunner()

    def _main_command_loop(self):
        create_query_db = CreateDBCommand(options={
            "input_files": self.query,
            "db_out_path": "db/query",
        })
        create_target_db = CreateDBCommand(options={
            "input_files": self.target,
            "db_out_path": "db/target",
        })

        search = SearchCommand(options={
            "query": "db/query",
            "target": "db/target",
            "result": "db/result",
            "search_tmp": "db/search_tmp",
            "a": "1",
        })

        convert_alis = ConvertAlisCommand(options={
            "query": "db/query",
            "target": "db/target",
            "result": "db/result",
            "output": "x_aln_x_.8",
            "format_mode": "0",
        })

        self.runner.add_command(create_query_db)
        self.runner.add_command(create_target_db, dependencies=[create_query_db])
        self.runner.add_command(search, dependencies=[create_target_db])
        self.runner.add_command(convert_alis, dependencies=[search])

    def seek(self):
        self._main_command_loop()
        self.runner.run()

    @property
    def output(self):
        if self.runner.is_successful():
            resuls = ""
            for output in self.runner.get_state().values():
                resuls += output + "\n"
            return resuls

        else:
            return f"Workflow failed with error: \n{self.runner.get_error()}"


class SimpleBashCommandsWorkflow:
    def __init__(self):
        self.runner = WorkflowRunner()

    def _main_command_loop(self):
        self.runner.add_command(
            BaseCommand("df", [""], {"help": ""})
        )

    def seek(self):
        self._main_command_loop()
        self.runner.run()

    @property
    def output(self):
        if self.runner.is_successful():
            resuls = ""
            for output in self.runner.get_state().values():
                resuls += output + "\n"
            return resuls

        else:
            return f"Workflow failed with error: \n{self.runner.get_error()}"

if __name__ == "__main__":
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

    w = SimpleBashCommandsWorkflow()
    w.seek()
    print(w.output)

    # runner.add_command(FoldSeekCommand("easy-search", [""], {"h": ""}))
