"""Pytest configuration file."""
from pathlib import Path

import pytest
from _pytest.config import Config
from _pytest.config.argparsing import Parser
from _pytest.fixtures import FixtureRequest


def pytest_configure(config: Config) -> None:
    """Add custom markers to pytest.

    Args:
        config (Config): The pytest config.
    """
    config.addinivalue_line("markers", "contacts: mark test as contacts test")
    config.addinivalue_line("markers", "ag: mark test as MDAnalysis AtomGroup support test")
    config.addinivalue_line("markers", "nb: mark test as NeighborPairs test")
    config.addinivalue_line("markers", "au: mark test as array_utils test")
    config.addinivalue_line("markers", "trajs: mark test as trajectory test")


def pytest_addoption(parser: Parser) -> None:
    """Add command line options to pytest.

    Args:
        parser (Parser): The pytest parser.
    """
    parser.addoption("--large-files", action="store_true", default=False, help="Use large files for testing")
    parser.addoption("--contacts", action="store_true", default=False, help="Run contacts tests")
    parser.addoption("--ag", action="store_true", default=False, help="Run MDAnalysis AtomGroup support tests")
    parser.addoption("--nb", action="store_true", default=False, help="Run NeighborPairs tests")
    parser.addoption("--au", action="store_true", default=False, help="Run array_utils tests")
    parser.addoption("--trajs", action="store_true", default=False, help="Run trajectory tests")


def pytest_collection_modifyitems(config: Config, items: list[pytest.Item]) -> None:
    """Modify the list of tests to run based on command line options.

    Args:
        config (Config): The pytest config.
        items (list[pytest.Item]): The list of tests to run.
    """
    if config.getoption("--contacts"):
        items[:] = [item for item in items if item.get_closest_marker("contacts")]
    elif config.getoption("--ag"):
        items[:] = [item for item in items if item.get_closest_marker("ag")]
    elif config.getoption("--nb"):
        items[:] = [item for item in items if item.get_closest_marker("nb")]
    elif config.getoption("--au"):
        items[:] = [item for item in items if item.get_closest_marker("au")]
    elif config.getoption("--trajs"):
        items[:] = [item for item in items if item.get_closest_marker("trajs")]


@pytest.fixture(scope="session", autouse=True)
def cleanup() -> None:  # noqa: PT004
    """Remove hidden files after testing. Specifically, MDAnalysis creates hidden files that end with .lock or .npz
    in the directory where the tests are run. This fixture removes those files after testing.

    Args:
        request (FixtureRequest): The pytest fixture request.
    """

    def remove_hidden_files() -> None:
        test_dir = Path(__file__).parent / "tests" / "data"

        # Find hidden files that end with .lock or .npz in that directory
        for pattern in (".*.lock", ".*.npz"):
            hidden_files = test_dir.glob(pattern)
            for file_path in hidden_files:
                if file_path.is_file():
                    file_path.unlink()

    def remove_downloaded_files() -> None:
        test_dir = Path(__file__).parent

        # Find downloaded files in that directory
        for pattern in ("*.pdb", "*.pdb.gz", "*.cif", "*.cif.gz"):
            downloaded_files = test_dir.glob(pattern)
            for file_path in downloaded_files:
                if file_path.is_file():
                    file_path.unlink()

    yield
    remove_hidden_files()
    remove_downloaded_files()
