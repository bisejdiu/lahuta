"""GROMACS is cool. Let's try to be cool, too."""
import re
from pathlib import Path


def parse_quotes(file_path: Path) -> list[tuple[str, str]]:
    """Parse quotes from the GROMACS coolstuff.cpp file.

    Args:
        file_path (str): Path to the coolstuff.cpp file.

    Returns:
        list[tuple[str, str]]: List of tuples containing the quote and the author.

    """
    with Path(file_path).open("r") as file:
        content = file.read()
        matches = re.findall(r'{ "(.*?)", "(.*?)" }', content)
        return [(match[0], match[1]) for match in matches]


PROJECT_ROOT = Path(__file__).parent.parent.parent / "tests" / "data"
GROMACS_QUOTES = parse_quotes(PROJECT_ROOT / "coolstuff.cpp")
