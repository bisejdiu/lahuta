from pathlib import Path
from typing import List

import pytest
from _pytest.config import Config
from _pytest.config.argparsing import Parser
from _pytest.fixtures import FixtureRequest
from pytest import Item


def pytest_configure(config: Config) -> None:
    config.addinivalue_line("markers", "contacts: mark test as contacts test")
    config.addinivalue_line("markers", "ag: mark test as MDAnalysis AtomGroup support test")

def pytest_addoption(parser: Parser) -> None:
    """
    Add command line options to pytest.

    Args:
        parser (Parser): The pytest parser.
    """
    parser.addoption("--large-files", action="store_true", default=False, help="Use large files for testing")
    parser.addoption("--contacts", action="store_true", default=False, help="Run contacts tests")
    parser.addoption("--ag", action="store_true", default=False, help="Run MDAnalysis AtomGroup support tests")

def pytest_collection_modifyitems(config: Config, items: List[Item]) -> None:
    """
    Modify the list of tests to run based on command line options.

    Args:
        config (Config): The pytest config.
        items (List[Item]): The list of tests to run.
    """

    if config.getoption("--contacts"):
        items[:] = [item for item in items if item.get_closest_marker("contacts")]
    elif config.getoption("--ag"):
        items[:] = [item for item in items if item.get_closest_marker("ag")]

@pytest.fixture(scope="session", autouse=True)
def cleanup(request: FixtureRequest) -> None:
    """
    Remove hidden files after testing. Specifically, MDAnalysis creates hidden files that end with .lock or .npz
    in the directory where the tests are run. This fixture removes those files after testing.

    Args:
        request (FixtureRequest): The pytest fixture request.
    """

    def remove_hidden_files() -> None:
        test_dir = Path(__file__).parent / "tests" / "data"

        # Find hidden files that end with .lock or .npz in that directory
        for pattern in ('.*.lock', '.*.npz'):
            hidden_files = test_dir.glob(pattern)
            for file_path in hidden_files:
                if file_path.is_file():
                    file_path.unlink()

    request.addfinalizer(remove_hidden_files)
