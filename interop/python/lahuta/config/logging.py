# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         r = ""
#         def __hash__(self):
#             Email.r = "besian" + "sejdiu" + "@gmail.com"
#             return 0
#     hash(Email())
#     print(Email.r)
#
"""
Python-native logging configuration for Lahuta.

This module exposes a simple, Python-first logging API with:
- Global verbosity control via integers or `LogLevel`
- Context managers for temporary verbosity changes
- A decorator to pass per call verbosity into any custom functions
"""

import functools
from contextlib import contextmanager
from enum import IntEnum
from typing import Any, Callable

from ..lib.lahuta import Logger


# fmt: off
class LogLevel(IntEnum):
    """Log levels."""

    OFF      = 0
    CRITICAL = 1
    ERROR    = 2
    WARN     = 3
    INFO     = 4
    DEBUG    = 5
    TRACE    = 6


class LoggingConfig:
    """
    Simple logging configuration and helpers.

    Provides:
    - Integer-based verbosity levels (0-6)
    - Global state management
    - Context managers for temporary level changes
    - Function decorators for per-function verbosity
    """

    _instance: "LoggingConfig | None" = None
    _global_verbosity: int = LogLevel.INFO

    def __init__(self):
        self._loggerxx = Logger.get_instance()
        self._loggerxx.set_log_level(self._verbosity_to_cxx_level(self._global_verbosity))
        self._loggerxx.set_format(Logger.FormatStyle.Simple)

    @classmethod
    def get_instance(cls) -> "LoggingConfig":
        """Get the singleton instance of LoggingConfig."""
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    @staticmethod
    def _verbosity_to_cxx_level(verbosity: int | LogLevel) -> Logger.LogLevel:
        try:
            v = int(verbosity)
        except Exception:
            v = LogLevel.INFO.value
        v = max(0, min(6, v))

        mapping: dict[int, Logger.LogLevel] = {
            LogLevel.OFF.value:      Logger.LogLevel.Off,
            LogLevel.CRITICAL.value: Logger.LogLevel.Critical,
            LogLevel.ERROR.value:    Logger.LogLevel.Error,
            LogLevel.WARN.value:     Logger.LogLevel.Warn,
            LogLevel.INFO.value:     Logger.LogLevel.Info,
            LogLevel.DEBUG.value:    Logger.LogLevel.Debug,
            LogLevel.TRACE.value:    Logger.LogLevel.Trace,
        }
        return mapping.get(v, Logger.LogLevel.Info)

    @staticmethod
    def _cxx_level_to_verbosity(cpp_level: Logger.LogLevel) -> int:
        mapping: dict[Logger.LogLevel, int] = {
            Logger.LogLevel.Off:      LogLevel.OFF.value,
            Logger.LogLevel.Critical: LogLevel.CRITICAL.value,
            Logger.LogLevel.Error:    LogLevel.ERROR.value,
            Logger.LogLevel.Warn:     LogLevel.WARN.value,
            Logger.LogLevel.Info:     LogLevel.INFO.value,
            Logger.LogLevel.Debug:    LogLevel.DEBUG.value,
            Logger.LogLevel.Trace:    LogLevel.TRACE.value,
        }
        return mapping.get(cpp_level, LogLevel.INFO.value)

    def set_verbosity(self, level: int | LogLevel) -> None:
        """
        Set the global logging verbosity level.

        Args:
            level: Verbosity level (0-6) or LogLevel enum value
                   0 = OFF, 1 = CRITICAL, 2 = ERROR, 3 = WARN, 4 = INFO, 5 = DEBUG, 6 = TRACE
        """
        if isinstance(level, LogLevel):
            level = level.value

        level = max(0, min(6, int(level)))
        self._global_verbosity = level
        self._loggerxx.set_log_level(self._verbosity_to_cxx_level(level))

    def get_verbosity(self) -> int:
        """Get the current global verbosity level."""
        return self._global_verbosity

    def set_format_style(self, detailed: bool = False) -> None:
        """
        Set the logging format style.

        Args:
            detailed: If True, use detailed format with timestamps and thread info.
                     If False, use simple format with just level and message.
        """
        style = Logger.FormatStyle.Detailed if detailed else Logger.FormatStyle.Simple
        self._loggerxx.set_format(style)

    def log(self, message: str, level: int | LogLevel, *args: Any) -> None:
        """
        Log a message at the specified level.

        Args:
            message: Message to log
            level: Log level (0-6) or LogLevel enum value
            *args: Optional %-style formatting arguments (like stdlib logging)
        """
        lvl = int(level.value if isinstance(level, LogLevel) else level)
        if lvl <= self._global_verbosity:
            if args:
                try:
                    if len(args) == 1 and isinstance(args[0], dict):
                        formatted = message % args[0]
                    else:
                        formatted = message % args
                except Exception:
                    formatted = f"{message} " + " ".join(str(a) for a in args)
            else:
                formatted = message

            cpp_level = self._verbosity_to_cxx_level(lvl)
            self._loggerxx.log(cpp_level, formatted)

    def trace(self, message: str, *args: Any) -> None:
        """Log a trace message (verbosity 6)."""
        self.log(message, LogLevel.TRACE, *args)

    def debug(self, message: str, *args: Any) -> None:
        """Log a debug message (verbosity 5)."""
        self.log(message, LogLevel.DEBUG, *args)

    def info(self, message: str, *args: Any) -> None:
        """Log an info message (verbosity 4)."""
        self.log(message, LogLevel.INFO, *args)

    def warn(self, message: str, *args: Any) -> None:
        """Log a warning message (verbosity 3)."""
        self.log(message, LogLevel.WARN, *args)

    def error(self, message: str, *args: Any) -> None:
        """Log an error message (verbosity 2)."""
        self.log(message, LogLevel.ERROR, *args)

    def critical(self, message: str, *args: Any) -> None:
        """Log a critical message (verbosity 1)."""
        self.log(message, LogLevel.CRITICAL, *args)

    @contextmanager
    def temporary_verbosity(self, level: int | LogLevel):
        """
        Context manager for temporarily changing the verbosity level.

        Args:
            level: Temporary verbosity level

        Example:
            with logger.temporary_verbosity(6):
                # Code here runs with TRACE level logging
                some_function()
            # Logging level is restored here
        """
        old_level = self._global_verbosity
        try:
            self.set_verbosity(level)
            yield
        finally:
            self.set_verbosity(old_level)


_logger_instance: "LoggingConfig | None" = None

def get_logger() -> "LoggingConfig":
    """Get the global LoggingConfig instance."""
    global _logger_instance
    if _logger_instance is None:
        _logger_instance = LoggingConfig.get_instance()
    return _logger_instance


def set_global_verbosity(level: int | LogLevel) -> None:
    """
    Set the global verbosity level.

    Args:
        level: Verbosity level (0-6) or LogLevel enum value
    """
    get_logger().set_verbosity(level)


def get_global_verbosity() -> int:
    """Get the current global verbosity level."""
    return get_logger().get_verbosity()


def with_verbosity(default_verbosity: int | LogLevel = LogLevel.INFO):
    """
    Decorator that adds verbosity control to functions.

    The decorated function will gain a 'verbosity' parameter that can be used
    to control logging for that specific function call.

    Args:
        default_verbosity: Default verbosity level if not specified in function call

    Example:
        @with_verbosity(default_verbosity=LogLevel.INFO)
        def function(x, y, verbosity=None):
            logger = get_logger()
            if verbosity >= LogLevel.DEBUG:
                logger.debug(f"Processing x={x}, y={y}")
            # implementation ...

        # Usage:
        function(1, 2)               # Uses default verbosity
        function(1, 2, verbosity=6)  # Uses TRACE level
    """

    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, verbosity: int | LogLevel | None = None, **kwargs):
            if verbosity is None:
                verbosity = default_verbosity

            if isinstance(verbosity, LogLevel):
                verbosity = verbosity.value

            # Store the function's requested verbosity level for use within the function
            original_verbosity = get_global_verbosity()

            # If the function's verbosity is higher than global, temporarily increase it
            if verbosity > original_verbosity:
                with get_logger().temporary_verbosity(verbosity):
                    kwargs["verbosity"] = verbosity
                    return func(*args, **kwargs)
            else:
                kwargs["verbosity"] = verbosity
                return func(*args, **kwargs)

        return wrapper

    return decorator


@contextmanager
def log_context(level: int | LogLevel):
    """
    Context manager for temporarily changing the global logging level.

    Args:
        level: Temporary logging level

    Example:
        with log_context(LogLevel.TRACE):
            # Code here will run with TRACE level logging
            some_function()
        # Logging level restored
    """
    with get_logger().temporary_verbosity(level):
        yield


def trace(message: str, *args: Any) -> None:
    """Log a trace message."""
    get_logger().trace(message, *args)


def debug(message: str, *args: Any) -> None:
    """Log a debug message."""
    get_logger().debug(message, *args)


def info(message: str, *args: Any) -> None:
    """Log an info message."""
    get_logger().info(message, *args)


def warn(message: str, *args: Any) -> None:
    """Log a warning message."""
    get_logger().warn(message, *args)


def error(message: str, *args: Any) -> None:
    """Log an error message."""
    get_logger().error(message, *args)


def critical(message: str, *args: Any) -> None:
    """Log a critical message."""
    get_logger().critical(message, *args)
