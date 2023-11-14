import asyncio
from dep_graph import DependencyGraph
from cmd_base import FoldSeekCommand
from createdb import CreateDBCommand

class WorkflowRunner:
    def __init__(self):
        self.dependency_graph = DependencyGraph()
        self.state = {}
        self.success = True
        self.error = None

    def add_command(self, command, dependencies=[]):
        self.dependency_graph.add_command(command, dependencies=dependencies)

    def run(self):
        try:
            loop = asyncio.get_event_loop()
            self.state = loop.run_until_complete(self.dependency_graph.execute_graph())
        except Exception as e:
            self.success = False
            self.error = str(e)

    def is_successful(self):
        return self.success

    def get_state(self):
        return self.state

    def get_error(self):
        return self.error


runner = WorkflowRunner()

# createdb = CreateDBCommand(options={
#     "input_files": "1gzm.pdb",
#     "db_out_path": "db/query",
#     "chain_name_mode": "0",
# })
# createdb2 = CreateDBCommand(options={
#     "input_files": "examples/",
#     "db_out_path": "db/target",
#     "chain_name_mode": "0",
# })

runner.add_command(FoldSeekCommand("easy-search", [""], {"h": ""}))
# runner.add_command(createdb)
# runner.add_command(createdb2, dependencies=[createdb])

runner.run()

if runner.is_successful():
    print("Workflow completed successfully.")
    for output in runner.get_state().values():
        print(output)
else:
    print("Workflow failed with error:", runner.get_error())