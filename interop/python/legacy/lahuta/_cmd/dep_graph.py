from typing import Optional

from .base import BaseCommand


class DependencyGraph:
    def __init__(self) -> None:
        self.graph: dict[BaseCommand, list[BaseCommand]] = {}  # Stores the graph
        self.visited: dict[BaseCommand, bool] = {}  # Tracks visited nodes
        self.stack: dict[BaseCommand, bool] = {}  # Tracks nodes in the current DFS path

    def add_command(self, command: BaseCommand, dependencies: Optional[list[BaseCommand]] = None) -> None:
        self.graph[command] = dependencies or []
        self.visited[command] = False
        self.stack[command] = False

    async def execute_command(self, command: BaseCommand) -> str:
        # Run the command and handle the result
        success = await command.run()
        if not success:
            raise RuntimeError(f"Command {command} failed.", command.error)

        return command.output

    async def execute_graph(self) -> dict[str, str]:
        # Topologically sort the graph and execute commands
        results = {}
        for command in self.topological_sort():
            output = await self.execute_command(command)
            results[f"{command}"] = output
        return results

    def topological_sort(self) -> list[BaseCommand]:
        order: list[BaseCommand] = []
        for command in self.graph:
            if not self.visited[command] and self._topological_sort_util(command, order):
                raise Exception("Cycle detected in dependencies")  # noqa: TRY002
        return order

    def _topological_sort_util(self, command: BaseCommand, order: list[BaseCommand]) -> bool:
        self.visited[command] = True
        self.stack[command] = True

        for dep in self.graph[command]:
            if not self.visited[dep]:
                if self._topological_sort_util(dep, order):
                    return True
            elif self.stack[dep]:
                return True  # Cycle detected

        self.stack[command] = False
        order.append(command)
        return False
