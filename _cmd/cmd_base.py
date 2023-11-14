import asyncio

class Command:
    def __init__(self, command, args, kwargs):
        self.command = command
        # self._a = {f"--{key.replace('_', '-')} {value}" for key, value in kwargs.items()}
        self.args = []
        for key, value in kwargs.items():
            arg = f"-{key}" if len(key) == 1 else f"--{key.replace('_', '-')}"
            self.args.extend([arg, str(value)])

        self.args = tuple(args) + tuple(self.args)
        self.output = None
        self.error = None
        # print("type: ", self.args, type(self.args)) 

    async def run(self):
        print(f"Running command: {self}")
        process = await asyncio.create_subprocess_exec(
            "foldseek", self.command, *self.args,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE)
        
        stdout, stderr = await process.communicate()
        if process.returncode == 0:
            self.output = stdout.decode().strip()
        else:
            self.error = stderr.decode().strip()

        return process.returncode == 0
    
    @property
    def name(self):
        return self.command
    
    @property
    def arguments(self):
        return self.args

    def __str__(self):
        return f"{self.command} {' '.join(self.args)}"

