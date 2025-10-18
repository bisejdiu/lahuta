"""Demonstrates the logging capabilities of the lahuta library."""

from lahuta import logging
from lahuta.logging import LogLevel


def logging_quickstart() -> None:
    logging.set_global_verbosity(LogLevel.INFO)
    logging.info("hello from lahuta logging")
    logging.info("percent style: %s" % "works")
    logging.debug("this debug is hidden at INFO level")
    logging.warn("This is a warning message")

    with logging.log_context(LogLevel.DEBUG):
        logging.debug("debug visible inside context")
        logging.trace("trace still hidden (level=DEBUG)")

    logging.info("done")


def logging_with_decorator() -> int:
    """Demonstrates the use of a logging decorator to set function verbosity."""

    @logging.with_verbosity(default_verbosity=LogLevel.INFO)
    def compute(n: int, verbosity=None) -> int:
        if verbosity is None:
            verbosity = logging.get_global_verbosity()
        total = 0
        for i in range(n):
            if verbosity >= LogLevel.DEBUG and i % 10 == 0:
                logging.debug(f"i={i} total={total}")
            total += i * i
        return total

    logging.set_global_verbosity(LogLevel.WARN)
    a = compute(25)  # uses default INFO internally
    b = compute(25, verbosity=LogLevel.DEBUG)
    logging.info(f"compute results: {a}, {b}")
    return b


if __name__ == "__main__":
    logging_quickstart()
    result = logging_with_decorator()
    logging.info(f"Final result from compute: {result}")
