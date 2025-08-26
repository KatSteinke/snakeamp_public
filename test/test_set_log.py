import logging
import pathlib
import re
import unittest

import pytest

import set_log


class TestRecordList(unittest.TestCase):
    def test_record_list(self):
        """Save log entries individually."""
        test_handler_1 = set_log.RecordsListHandler()
        test_logger = logging.getLogger('test_record_list')
        test_logger.setLevel(logging.INFO)
        test_logger.addHandler(test_handler_1)
        assert not test_handler_1.record_log
        test_logger.info('Testing 1')
        assert len(test_handler_1.record_log) == 1
        assert test_handler_1.record_log[0].name == "test_record_list"
        assert test_handler_1.record_log[0].levelno == logging.INFO
        assert test_handler_1.record_log[0].message == "Testing 1"
        test_handler_2 = set_log.RecordsListHandler()
        test_logger.addHandler(test_handler_2)
        assert not test_handler_2.record_log
        test_logger.info('Testing 2')
        assert len(test_handler_2.record_log) == 1
        assert len(test_handler_1.record_log) == 2
        assert test_handler_2.record_log[0].name == "test_record_list"
        assert test_handler_2.record_log[0].levelno == logging.INFO
        assert test_handler_2.record_log[0].message == "Testing 2"


class TestHandleBulkLogs(unittest.TestCase):
    test_file = pathlib.Path(__file__).parent / "data" / "set_log" / "test.log"

    def setUp(self):
        with open(self.test_file, "w") as write_test:
            write_test.write("")

    def tearDown(self):
        with open(self.test_file, "w") as write_test:
            write_test.write("")

    def test_success_file_handler(self):
        """Write specified records to a file."""
        logs = [logging.LogRecord(name = "testlog", level = logging.INFO, pathname = "test/file",
                                  msg = "Testing 1", lineno = 1, args = {}, exc_info = None),
                logging.LogRecord(name = "testlog", level = logging.INFO, pathname = "test/file",
                                  msg = "Testing 2", lineno = 2, args = {}, exc_info = None)
                ]
        file_handler = logging.FileHandler(self.test_file)
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(logging.Formatter("%(name)s - %(levelno)s - %(message)s"))
        set_log.handle_bulk_logs(file_handler, logs)
        with open(self.test_file, "r") as read_test:
            lines = read_test.readlines()
            assert lines[0] == "testlog - 20 - Testing 1\n"
            assert lines[1] == "testlog - 20 - Testing 2\n"


class TestGetStreamLog(unittest.TestCase):
    def test_fail_bad_level(self):
        """Fail if an invalid log level was set."""
        error_msg = "Level FYI is not a valid log level."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            set_log.get_stream_log("bad_level", level="FYI")

    def test_success_default_level(self):
        """Get a stream logger with the default log level."""
        test_logger = set_log.get_stream_log("default")
        assert test_logger.name == "default"
        assert test_logger.level == 20  # INFO
        assert test_logger.hasHandlers()
        assert isinstance(test_logger.handlers[0], logging.StreamHandler)

    def test_success_set_level(self):
        """Get a stream logger with the specified log level."""
        test_logger = set_log.get_stream_log("debug_logger", level="DEBUG")
        assert test_logger.name == "debug_logger"
        assert test_logger.level == 10  # DEBUG
        assert test_logger.hasHandlers()
        assert isinstance(test_logger.handlers[0], logging.StreamHandler)

    def test_success_set_root(self):
        """Set stream logging on the root logger."""
        test_logger = set_log.get_stream_log()
        assert test_logger.name == "root"
        assert test_logger.level == 20
        assert test_logger.hasHandlers()
