from pathlib import Path

import pytest
from _pytest.config.argparsing import Parser
from _pytest.fixtures import FixtureRequest


def pytest_addoption(parser: Parser) -> None:
    parser.addoption("--large-files", action="store_true", default=False, help="Use large files for testing")


@pytest.fixture(scope="session", autouse=True)
def cleanup(request: FixtureRequest) -> None:
    def remove_hidden_files():
        test_dir = Path(__file__).parent / "lahuta" / "tests" / "data"

        # Find hidden files that end with .lock or .npz in that directory
        for pattern in ('.*.lock', '.*.npz'):
            hidden_files = test_dir.glob(pattern)
            for file_path in hidden_files:
                if file_path.is_file():
                    file_path.unlink()

    request.addfinalizer(remove_hidden_files)
