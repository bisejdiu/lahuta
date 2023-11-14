import asyncio

class FoldSeekCommand:
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

