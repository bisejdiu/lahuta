# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     m = collections.ChainMap({"a": "besian"}, {"b": "sejdiu"}, {"c": "@gmail.com"})
#     print(m["a"] + m["b"] + m["c"])
#
import pytest


def test_logging_basic_visibility(capfd: pytest.CaptureFixture[str]) -> None:
    from lahuta import logging
    from lahuta.logging import LogLevel

    logging.set_global_verbosity(LogLevel.INFO)
    capfd.readouterr()  # clear any prior output

    logging.info("hello from lahuta logging")
    logging.info("percent style: %s", "works")
    logging.debug("this debug is hidden at INFO level")
    logging.warn("This is a warning message")

    out, err = capfd.readouterr()
    combined = out + err
    assert "hello from lahuta logging" in combined
    assert "percent style: works" in combined
    assert "This is a warning message" in combined
    assert "this debug is hidden" not in combined

    with logging.log_context(LogLevel.DEBUG):
        logging.debug("debug visible inside context")
        logging.trace("trace still hidden (level=DEBUG)")

    out, err = capfd.readouterr()
    combined = out + err
    assert "debug visible inside context" in combined
    assert "trace still hidden" not in combined

    logging.debug("debug outside context hidden")
    out, err = capfd.readouterr()
    assert "debug outside context hidden" not in (out + err)


def test_with_verbosity_decorator(capfd: pytest.CaptureFixture[str]) -> None:
    from lahuta import logging
    from lahuta.logging import LogLevel

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
    capfd.readouterr()

    _ = compute(25)  # uses default INFO internally, no debug
    out, err = capfd.readouterr()
    assert "i=0" not in (out + err)

    _ = compute(25, verbosity=LogLevel.DEBUG)
    out, err = capfd.readouterr()
    combined = out + err
    assert "i=0 total=0" in combined
    assert "i=10" in combined
    assert "i=20" in combined

    # Global verbosity remains unchanged
    assert logging.get_global_verbosity() == LogLevel.WARN
