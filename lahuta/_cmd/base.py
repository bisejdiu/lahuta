import asyncio
from typing import Any, Mapping, Optional


async def run_command_async(command: str, sub_command: str, *args: str) -> tuple[str, str]:
    arguments = [command] + ([sub_command] if sub_command else []) + list(args)
    arguments = list(filter(lambda x: x != "", arguments))

    process = await asyncio.create_subprocess_exec(
        *arguments, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE
    )

    stdout, stderr = await process.communicate()
    if process.returncode == 0:
        return stdout.decode().strip(), ""

    return "", stderr.decode().strip()


async def run_command_async2(command: str, sub_command: str, *args: str) -> tuple[str, str]:
    arguments = [command] + ([sub_command] if sub_command else []) + list(args)
    arguments = list(filter(lambda x: x != "", arguments))

    process = await asyncio.create_subprocess_exec(
        *arguments, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE
    )

    if process.stdout is None or process.stderr is None:
        return "", ""

    stdout, stderr = [], []
    async for line in process.stdout:
        stdout.append(line.decode().strip())
    async for line in process.stderr:
        stderr.append(line.decode().strip())

    await process.wait()

    return "\n".join(stdout), "\n".join(stderr)


class BaseCommand:
    def __init__(
        self, command: str, sub_command: str, required_args: list[str], kwargs: Optional[Mapping[str, Any]] = None
    ) -> None:
        self.command = command
        self.sub_command = sub_command
        self.args = required_args + (self._handle_kwargs(kwargs) if kwargs else [])

        self.output = ""
        self.error = ""

    def _handle_kwargs(self, kwargs: Mapping[str, Any]) -> list[str]:
        arguments: list[str] = []
        for key, value in kwargs.items():
            arg = f"-{key}" if len(key) == 1 else f"--{key.replace('_', '-')}"
            arguments.extend([arg, str(value)])

        return arguments

    async def run(self) -> bool:
        raise NotImplementedError("BaseCommand.run() must be implemented by subclasses")

    @property
    def name(self) -> str:
        return self.sub_command

    @property
    def arguments(self) -> list[str]:
        return self.args

    def __str__(self) -> str:
        return f"{self.command} {self.sub_command} {' '.join(self.args)}"


class SimpleCommand(BaseCommand):
    def __init__(self, command: str, kwargs: Optional[Mapping[str, Any]] = None):
        super().__init__(command, "", [], kwargs)

    async def run(self) -> bool:
        self.output, self.error = await run_command_async2(self.command, self.sub_command, *self.args)
        return self.error == ""


class CommandWithoutSubCommand(BaseCommand):
    def __init__(self, command: str, required_args: list[str], kwargs: Optional[Mapping[str, Any]] = None):
        super().__init__(command, "", required_args, kwargs)

    async def run(self) -> bool:
        self.output, self.error = await run_command_async(self.command, self.sub_command, *self.args)
        return self.error == ""


class CommandWithSubCommand(BaseCommand):
    async def run(self) -> bool:
        self.output, self.error = await run_command_async(self.command, self.sub_command, *self.args)
        return self.error == ""


class TestBaseCommand(BaseCommand):
    NAME = "testcommand"

    def __init__(self, command: str, sub_command: str, required_args: list[str], kwargs: Mapping[str, Any]) -> None:
        super().__init__(command, sub_command, required_args, kwargs)
        self.returncode = 0

    async def run(self) -> bool:
        arguments = [self.NAME] + ([self.sub_command] if self.sub_command else []) + self.args
        arguments = list(filter(lambda x: x != "", arguments))

        print("Running command:", *arguments)

        if self.returncode == 0:
            self.output = self.output
        else:
            self.error = self.error

        return self.returncode == 0


async def main() -> None:
    test_cmd = TestBaseCommand("abc", "some_command", ["arg1", "arg2"], {"option": "value"})
    test_cmd.output = "Expected output for success"
    test_cmd.error = "Simulated error message"
    test_cmd.returncode = 0

    # Run the command and assert based on the output/error
    await test_cmd.run()
    assert test_cmd.output == "Expected output for success", test_cmd.output


if __name__ == "__main__":
    asyncio.run(main())
