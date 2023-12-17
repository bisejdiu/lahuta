"""Downloader."""
import asyncio
import logging
from enum import Enum
from pathlib import Path

import httpx
import rich.progress
from rich.progress import Progress, TaskID


class URL(Enum):
    """URLs."""

    PDB = "https://files.rcsb.org/download/"


class FileDownloader:
    """Class to download files from a URL.

    Given a list of file names, a URL, and a directory, it downloads the files from the URL and saves them in the
    directory. If the directory does not exist, it creates it (along with other sanity checks).

    Attributes:
        file_names (list[str]): List of file names.
        url (str): URL.
        dir_loc (Path): Path to the directory.

    Methods:
        download_all: Download all files.
    """

    def __init__(self, *, file_names: list[str], dir_loc: str | Path | None = None, url: str | URL = URL.PDB):
        self.file_names = file_names
        self.url = url.value if isinstance(url, URL) else url
        self.dir_loc = Path(dir_loc) if dir_loc else Path.cwd()
        self.dir_loc.mkdir(parents=True, exist_ok=True)

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

    async def _download_file(
        self, client: httpx.AsyncClient, file_name: str, progress: Progress, download_task: TaskID
    ) -> None:
        local_path = self.dir_loc / file_name
        if local_path.exists():
            progress.update(download_task, advance=1)
            return

        response = await client.get(f"{self.url}{file_name}", follow_redirects=True)
        if response.status_code == 200:
            with local_path.open("wb") as f:
                f.write(response.content)
            progress.update(download_task, advance=1)

    async def _download_all_async(self) -> None:
        limits = httpx.Limits(max_keepalive_connections=20, max_connections=10)
        transport = httpx.AsyncHTTPTransport(retries=1, http2=True)
        async with httpx.AsyncClient(transport=transport, limits=limits, timeout=None) as client:
            with rich.progress.Progress(
                "[progress.percentage]{task.percentage:>3.0f}%",
                rich.progress.BarColumn(bar_width=None),
                rich.progress.TextColumn("{task.completed}/{task.total} files"),
                rich.progress.TextColumn("[progress.elapsed]Elapsed time: {task.elapsed:.4f}"),
            ) as progress:
                download_task = progress.add_task("Downloading files...", total=len(self.file_names))
                await asyncio.gather(
                    *[self._download_file(client, name, progress, download_task) for name in self.file_names]
                )

    def download_all(self) -> None:
        """Download all files."""
        asyncio.run(self._download_all_async())


def easy_downloader(*, url: str | URL | None, file_names: list[str], dir_loc: str | Path | None = None) -> None:
    """Easy downloader.

    A convenient wrapper around the FileDownloader class. It takes a list of file names, a URL, and a directory, and
    downloads the files from the URL to the directory. If the directory does not exist, it creates it (along with other
    sanity checks). If the URL is None, it defaults to the PDB URL.

    Args:
        url (str | URL | None): URL.
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
        url = URL.PDB

    downloader = FileDownloader(file_names=file_names, dir_loc=dir_loc, url=url)
    downloader.download_all()


# Usage
if __name__ == "__main__":
    cif_path = "/home/bisejdiu/projects/landscapy/data/downloaded/cifs"

    cif_codes = {}
    for cif in Path(cif_path).glob("*.cif"):
        cif_codes[cif.name] = cif

    sample_cif_codes = list(cif_codes.keys())  # [:100]
    print(len(sample_cif_codes))

    easy_downloader(url=URL.PDB, file_names=sample_cif_codes, dir_loc=Path("test1/test2"))
