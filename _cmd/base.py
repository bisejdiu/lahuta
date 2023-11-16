import asyncio
from typing import Any, Iterable, Mapping

class BaseCommand:
    NAME = ""
    def __init__(self, sub_command: str, required_args: Iterable[str], kwargs: Mapping[str, Any]):
        self.sub_command = sub_command
        self.args = required_args + self._handle_kwargs(kwargs)

        self.output = ""
        self.error = ""

    def _handle_kwargs(self, kwargs):
        arguments = []
        for key, value in kwargs.items():
            arg = f"-{key}" if len(key) == 1 else f"--{key.replace('_', '-')}"
            arguments.extend([arg, str(value)])

        return arguments

    async def run(self):
        raise NotImplementedError("BaseCommand.run() must be implemented by subclasses")
    
    @property
    def name(self):
        return self.sub_command
    
    @property
    def arguments(self):
        return self.args

    def __str__(self):
        return f"{self.NAME} {self.sub_command} {' '.join(self.args)}"

class FoldSeekBaseCommand(BaseCommand):
    NAME = "foldseek"

    async def run(self):
        arguments = [self.NAME] + ([self.sub_command] if self.sub_command else []) + self.args
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
    
class TestBaseCommand(BaseCommand):
    NAME = "testcommand"
    def __init__(self, sub_command: str, required_args: Iterable[str], kwargs: Mapping[str, Any]) -> None:
        super().__init__(sub_command, required_args, kwargs)
        self.returncode = 0

    async def run(self) -> bool:
        arguments = [self.NAME] + ([self.sub_command] if self.sub_command else []) + self.args
        arguments = list(filter(lambda x: x != '', arguments))

        if self.returncode == 0:
            self.output = self.output
        else:
            self.error = self.error

        return self.returncode == 0

async def main():
    test_cmd = TestBaseCommand("some_command", ["arg1", "arg2"], {"option": "value"})
    test_cmd.output = "Expected output for success"
    test_cmd.error = "Simulated error message"
    test_cmd.returncode = 0 

    # Run the command and assert based on the output/error
    await test_cmd.run()
    assert test_cmd.output == "Expected output for success", test_cmd.output


if __name__ == '__main__':
    asyncio.run(main())
