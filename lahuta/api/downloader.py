"""Downloader."""
import asyncio
import logging

# from memory_profiler import profile
import os
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Coroutine, Generator, Iterator, TypeVar, cast

import httpx
import tqdm

# type aliases
NetworkError = str
FilePath = tuple[str, Path | None]
Result = tuple[str, Path]
NoSavedResult = tuple[None, None]
httpxResult = FilePath | NoSavedResult
CoroutineResult = Coroutine[Any, Any, httpxResult]
TaskType = Generator[CoroutineResult, None, None]


class URLs(Enum):
    """URLs."""

    RCSB = "https://files.rcsb.org/download/"
    PDBeKB = "https://www.ebi.ac.uk/pdbe/entry-files/download/"
    AlphaFold = "https://alphafold.ebi.ac.uk/files/"


class ProgressBarType(Enum):
    """The type of progress bar to show."""

    RICH = 1
    TQDM = 2


T = TypeVar("T", bound="NonNullResultIterator")


class NonNullResultIterator:
    """Iterator that skips over None results.

    Args:
        result (Iterator[tuple[str | None, Path | None]]): Iterator to wrap.
    """

    def __init__(self, result: Iterator[tuple[str | None, Path | None]]):
        self.result = result

    def __iter__(self: T) -> T:
        return self

    def __next__(self) -> tuple[str, Path]:
        while True:
            file_name, file_path = next(self.result)
            if file_name is not None and file_path is not None:
                return file_name, file_path
            if file_name is None and file_path is None:
                raise StopIteration


class FileDownloader:
    """Class to download files from a URL.

    Given a list of file names, a URL, and a directory, it downloads the files from the URL and saves them in the
    directory. If the directory does not exist, it creates it (along with other sanity checks).

    Attributes:
        file_names (list[str]): List of file names.
        url (str): URL.
        dir_loc (Path): Path to the directory.
        progress_bar_type (ProgressBarType): Type of progress bar to show.
        list_generator_limit (int): Limit for the number of files before we require a generator.
        semaphore_limit (int): Limit for the number of concurrent downloads.

    Methods:
        download_all: Download all files.
    """

    def __init__(  # noqa: PLR0913
        self,
        *,
        file_names: list[str] | Generator[str, None, None],
        dir_loc: str | Path | None = None,
        url: str | URLs = URLs.RCSB,
        list_generator_limit: int = 10_000,
        semaphore_limit: int = 100,
        progress_bar_type: ProgressBarType = ProgressBarType.RICH,
    ):
        if not isinstance(file_names, Generator) and len(file_names) > list_generator_limit:
            raise ValueError(f"File number is larger than {list_generator_limit}, please use a generator.")

        self.total_length = len(file_names) if isinstance(file_names, list) else None
        self.file_names = (x for x in file_names)

        self.url = url.value if isinstance(url, URLs) else url
        self.dir_loc = Path(dir_loc) if dir_loc else Path.cwd()
        self.dir_loc.mkdir(parents=True, exist_ok=True)
        self.progress_bar_type = progress_bar_type
        self.semaphore_limit = semaphore_limit

        self._check_health()

    def _check_health(self) -> None:
        if not self.file_names:
            raise ValueError("No file names provided.")

        if not self.url:
            raise ValueError("No URL provided.")

        if not self.dir_loc:
            raise ValueError("No directory provided.")

        if not self.dir_loc.exists():
            raise ValueError("Directory does not exist.")

        if not self.dir_loc.is_dir():
            raise ValueError("Directory is not a directory.")

        if any(self.dir_loc.iterdir()):
            logging.warning("Directory is not empty.")

        try:
            test_file_name = "lahuta_" + self.__class__.__name__ + "_test"
            test_file = self.dir_loc / test_file_name
            test_file.touch()
            test_file.unlink()
        except PermissionError:
            raise ValueError("Directory is not writable.") from None

    async def _generate_tasks(
        self, client: httpx.AsyncClient, progress_updater: Callable[[int], bool | None]
    ) -> TaskType:
        semaphore = asyncio.Semaphore(self.semaphore_limit)

        async def task_wrapper(name: str) -> httpxResult:
            async with semaphore:
                return await self._download_file(client, name, progress_updater)

        tasks = (task_wrapper(name) for name in self.file_names)
        return tasks

    async def _download_file(
        self, client: httpx.AsyncClient, file_name: str, progress_updater: Callable[[int], bool | None]
    ) -> httpxResult:
        local_path = self.dir_loc / file_name
        if local_path.exists():
            progress_updater(1)
            return file_name, local_path

        response = await client.get(f"{self.url}{file_name}", follow_redirects=True)
        if response.status_code == 200:
            with local_path.open("wb") as f:
                f.write(response.content)
            progress_updater(1)
            return file_name, local_path

        return file_name, None

    async def _download_files(self) -> tuple[list[Result], list[NetworkError]] | NoSavedResult:
        logging.getLogger("httpx").setLevel(logging.WARNING)
        limits = httpx.Limits(max_keepalive_connections=20, max_connections=10)
        transport = httpx.AsyncHTTPTransport(retries=1, http2=True)
        async with httpx.AsyncClient(transport=transport, limits=limits, timeout=None) as client:
            if self.progress_bar_type == ProgressBarType.RICH:
                result = await self._download_with_rich_progress(client)
            else:
                result = await self._download_with_tqdm_progress(client)

        output: list[Result] = []
        errors: list[NetworkError] = []
        for file_name, file_path in NonNullResultIterator(iter(result)):
            if file_path is None:
                errors.append(file_name)
                logging.warning(f"Could not download {file_name}.")
            else:
                output.append((file_name, file_path))

        return output, errors

    async def _download_with_rich_progress(self, client: httpx.AsyncClient) -> list[httpxResult]:
        import rich.progress

        def get_progress_bar_params(no_total: bool = False) -> tuple[Any, ...]:
            if no_total:
                # Style if no total value is known
                return (
                    "Download progress: ",
                    rich.progress.TextColumn("[bold blue]{task.completed}/? files", justify="right"),
                    rich.progress.SpinnerColumn(spinner_name="dots12", style="bold green", speed=3),
                    rich.progress.SpinnerColumn(spinner_name="pong", style="bold white", speed=1),
                )
            # Style parameters
            return (
                "[progress.percentage]{task.percentage:>3.0f}%",
                rich.progress.BarColumn(bar_width=None),
                rich.progress.TextColumn("{task.completed}/{task.total} files"),
                rich.progress.TextColumn("[progress.elapsed]Elapsed time: {task.elapsed:.4f}"),
                rich.progress.TimeRemainingColumn(),
            )

        bar_params = get_progress_bar_params(no_total=self.total_length is None)
        with rich.progress.Progress(*bar_params) as progress:
            download_task = progress.add_task("Downloading files...", total=self.total_length)
            tasks = await self._generate_tasks(client, lambda n: progress.update(download_task, advance=n))
            if self.total_length is None:
                await asyncio.gather(*tasks)
                return cast(list[httpxResult], [(None, None)])

            result = await asyncio.gather(*tasks)
            return cast(list[httpxResult], result)

    # @profile
    async def _download_with_tqdm_progress(self, client: httpx.AsyncClient) -> list[httpxResult]:
        with tqdm.tqdm(total=self.total_length, desc="Downloading files", unit=" files") as pbar:
            tasks: TaskType = await self._generate_tasks(client, lambda n: pbar.update(n))
            if self.total_length is None:
                await asyncio.gather(*tasks)
                return cast(list[httpxResult], [(None, None)])

            result = await asyncio.gather(*tasks)
            return cast(list[httpxResult], result)

    def download_all(self) -> asyncio.Task[tuple[list[Result], list[NetworkError]] | NoSavedResult]:
        """Download all files."""
        loop = asyncio.get_event_loop()
        task = loop.create_task(self._download_files())
        if not loop.is_running():
            loop.run_until_complete(task)
        return task


def fast_download(
    *, url: str | URLs | None, file_names: list[str] | Generator[str, None, None], dir_loc: str | Path | None = None
) -> asyncio.Task[tuple[list[Result], list[NetworkError]] | NoSavedResult]:
    """Easy downloader.

    A convenient wrapper around the FileDownloader class. It takes a list of file names, a URL, and a directory, and
    downloads the files from the URL to the directory. If the directory does not exist, it creates it (along with other
    sanity checks). If the URL is None, it defaults to the RCSB PDB URL.

    Args:
        url (str | URLs | None): URLs.
        file_names (list[str]): List of file names.
        dir_loc (str | Path | None): Path to the directory.
    """

    def is_valid_url(url: str) -> bool:
        """Check if the URL is valid."""
        from urllib.parse import urlparse

        parsed = urlparse(url)
        return bool(parsed.scheme) and bool(parsed.netloc)

    if isinstance(url, str) and not is_valid_url(url):
        raise ValueError("Invalid URL.")

    if url is None:
        url = URLs.RCSB

    downloader = FileDownloader(
        file_names=file_names,
        dir_loc=dir_loc,
        url=url,
        progress_bar_type=ProgressBarType.RICH,
        list_generator_limit=100_000,
    )
    task = downloader.download_all()
    return task


# Usage
if __name__ == "__main__":

    def file_generator(directory: str = "/mnt/f/PDB_ARCHIVE/mmCIF") -> Generator[str, None, None]:
        """Generate file names from a directory."""
        for _, _, files in os.walk(directory):
            for file in files:
                yield file

    cif_path = "/home/bisejdiu/projects/landscapy/data/downloaded/cifs"

    cif_codes = {}
    for cif in Path(cif_path).glob("*.cif"):
        cif_codes[cif.name] = cif

    sample_cif_codes = list(cif_codes.keys())  # [:500]

    fast_download(url=URLs.RCSB, file_names=file_generator(), dir_loc=Path("test"))
