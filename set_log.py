"""Convenience functions for logging."""

__author__ = "Kat Steinke"

import logging

from typing import List, Optional


class RecordsListHandler(logging.Handler):
    """A logging handler that stores log records for later use.
    Adapted from https://stackoverflow.com/a/62241832/

    Attributes:
        record_log:    the log records logged with this handler attached
    """
    def __init__(self):
        """Initialize the handler with a list for storing records."""
        super().__init__()
        self.record_log = []

    def emit(self, record):
        """Save a log record to the list.

        Arguments:
            record: the current log record

        """
        self.record_log.append(record)


def handle_bulk_logs(handler: logging.Handler, log_records: List[logging.LogRecord]) -> None:
    """Handle a number of log records with a given handler (i.e., save to file for a file handler,
    print for a stream handler,...)

    Arguments:
        handler:        the handler used for emitting the log records
        log_records:    the records to handle

    """
    for log_record in log_records:
        handler.emit(log_record)


def get_stream_log(logger_name: Optional[str]=None, level: str = "INFO") -> logging.Logger:
    """Create a stream logger with the specified level.

    Arguments:
        logger_name:    the name of the logger to set up
        level:          the logging level

    Returns:
        A logger with a StreamHandler attached, logging at the given level.

    Raises:
        ValueError: if a nonexistent logging level was given.
    """
    # 3.11 onwards: logging.getLevelNamesMapping()
    # - we'll hardcode for backwards compatibility for now
    log_levels = {"NOTSET": 0, "DEBUG": 10, "INFO": 20, "WARNING": 30, "ERROR": 40, "CRITICAL": 50}
    if level not in log_levels:
        raise ValueError(f"Level {level} is not a valid log level.")
    logger = logging.getLogger(logger_name)
    logger.setLevel(log_levels[level])
    console_log = logging.StreamHandler()
    console_log.setLevel(log_levels[level])
    logger.addHandler(console_log)
    return logger
