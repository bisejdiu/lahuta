import asyncio
from typing import Optional

from .base import BaseCommand
from .dep_graph import DependencyGraph


class WorkflowRunner:
    def __init__(self) -> None:
        self.dependency_graph = DependencyGraph()
        self.success = True
        self.error = ""

    def add_command(self, command: BaseCommand, dependencies: Optional[list[BaseCommand]] = None) -> None:
        dependencies = dependencies or []
        self.dependency_graph.add_command(command, dependencies=dependencies)

    async def run_async(self) -> dict[str, str]:
        try:
            return await self.dependency_graph.execute_graph()
        except Exception as e:  # noqa: BLE001
            self.success = False
            self.error = str(e)
            return {}

    def run(self) -> asyncio.Task[dict[str, str]]:
        loop = asyncio.get_event_loop()
        if loop.is_running():
            return asyncio.create_task(self.run_async())

        task = loop.create_task(self.run_async())
        loop.run_until_complete(task)
        return task

    def is_successful(self) -> bool:
        return self.success
