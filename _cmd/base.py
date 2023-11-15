import asyncio
from typing import Iterable

class BaseCommand:
    NAME = ""
    def __init__(self, command, required_args, kwargs):
        self.command = command
        self.args = required_args + self._handle_kwargs(kwargs)

        self.output = None
        self.error = None

    async def run(self):
        arguments = [self.NAME] + ([self.command] if self.command else []) + self.args
        arguments = list(filter(lambda x: x != '', arguments))

        process = await asyncio.create_subprocess_exec(
            *arguments,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE)
        
        stdout, stderr = await process.communicate()
        if process.returncode == 0:
            self.output = stdout.decode().strip()
        else:
            self.error = stderr.decode().strip()

        return process.returncode == 0
    
    def _handle_kwargs(self, kwargs):
        arguments = []
        for key, value in kwargs.items():
            arg = f"-{key}" if len(key) == 1 else f"--{key.replace('_', '-')}"
            arguments.extend([arg, str(value)])

        return arguments

    @property
    def name(self):
        return self.command
    
    @property
    def arguments(self):
        return self.args

    def __str__(self):
        return f"{self.NAME} {self.command} {' '.join(self.args)}"

class FoldSeekBaseCommand:
    NAME = "foldseek"
    def __init__(self, command, required_args, kwargs):
        self.command = command
        self.args = required_args + self._handle_kwargs(kwargs)

        self.output = None
        self.error = None

    async def run(self):
        arguments = [self.NAME] + ([self.command] if self.command else []) + self.args
        arguments = list(filter(lambda x: x != '', arguments))

        process = await asyncio.create_subprocess_exec(
            *arguments,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE)
        
        stdout, stderr = await process.communicate()
        if process.returncode == 0:
            self.output = stdout.decode().strip()
        else:
            self.error = stderr.decode().strip()

        return process.returncode == 0
    
    def _handle_kwargs(self, kwargs):
        arguments = []
        for key, value in kwargs.items():
            arg = f"-{key}" if len(key) == 1 else f"--{key.replace('_', '-')}"
            arguments.extend([arg, str(value)])

        return arguments

    @property
    def name(self):
        return self.command
    
    @property
    def arguments(self):
        return self.args

    def __str__(self):
        return f"{self.NAME} {self.command} {' '.join(self.args)}"

class TestBaseCommand:
    NAME = "testcommand"
    def __init__(self, command: str, required_args: Iterable[str], kwargs: dict[str, str]) -> None:
        self.command = command
        self.args = required_args + self._handle_kwargs(kwargs)

        # For testing, we can set these manually in tests
        self.output: str = ""
        self.error: str = ""
        self.returncode = 0

    async def run(self) -> bool:
        arguments = [self.NAME] + ([self.command] if self.command else []) + self.args
        arguments = list(filter(lambda x: x != '', arguments))

        if self.returncode == 0:
            self.output = self.output
        else:
            self.error = self.error

        return self.returncode == 0

    def _handle_kwargs(self, kwargs):
        arguments = []
        for key, value in kwargs.items():
            arg = f"-{key}" if len(key) == 1 else f"--{key.replace('_', '-')}"
            arguments.extend([arg, str(value)])

        return arguments

    @property
    def name(self):
        return self.command

    @property
    def arguments(self):
        return self.args

    def __str__(self):
        return f"{self.NAME} {self.command} {' '.join(self.args)}"

async def main():
    test_cmd = TestBaseCommand("some_command", ["arg1", "arg2"], {"option": "value"})
    test_cmd.output = "Expected output for success"
    test_cmd.error = "Simulated error message"
    test_cmd.returncode = 0  # Set to non-zero to simulate an error

    # Run the command and assert based on the output/error
    await test_cmd.run()
    assert test_cmd.output == "Expected output for success"


if __name__ == '__main__':
    asyncio.run(main())