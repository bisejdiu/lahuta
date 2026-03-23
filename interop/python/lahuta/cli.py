"""CLI wrapper that invokes the native lahuta binary."""

import os
import sys
from pathlib import Path


def main():
    """Execute the native lahuta CLI binary, replacing this process."""
    binary = Path(__file__).parent.parent / "bin" / "lahuta"
    os.execv(str(binary), [str(binary)] + sys.argv[1:])


if __name__ == "__main__":
    main()
