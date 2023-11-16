import asyncio
from dep_graph import DependencyGraph

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
