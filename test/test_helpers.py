import pathlib
import re
import unittest

import pytest

import helpers


class TestFindRundir(unittest.TestCase):
    def test_find_absolute_path_success(self):
        test_path = pathlib.Path(__file__).parent / "data" / "snake_helpers" \
                    / "test_dir_multi_pass" / "rawdata" / "subdir" / "fastq_pass"
        true_path = pathlib.Path(__file__).parent / "data" / "snake_helpers" \
                    / "test_dir_multi_pass" / "rawdata" / "subdir" / "fastq_pass"
        with self.assertLogs("helpers", level = "INFO") as logged:
            log_msg = f"INFO:helpers:Data is retrieved from the following folder:\n" \
                      f"{true_path}"
            test_fastq = helpers.get_fastq_pass_dir(test_path)
            assert log_msg in logged.output
        assert test_fastq == true_path

    def test_return_fastq_when_given(self):
        test_path = pathlib.Path(__file__).parent / "data" / "snake_helpers" / "test_dir"
        true_path = pathlib.Path(__file__).parent / "data" / "snake_helpers" / "test_dir" \
                    / "rawdata" / "subdir" / "fastq_pass"
        with self.assertLogs("helpers", level = "INFO") as logged:
            log_msg = f"INFO:helpers:Data is retrieved from the following folder:\n" \
                      f"{true_path}"
            test_fastq = helpers.get_fastq_pass_dir(test_path)
            assert log_msg in logged.output
        assert test_fastq == true_path

    def test_fail_path(self):
        test_path = pathlib.Path(__file__).parent / "data" / "helpers" / "subdir"
        error_msg = "fastq_pass folder(s) not found in expected location:\n" \
                    f"{test_path}/rawdata/*/fastq_pass\n" \
                    "Ensure correct directory and/or directory structure is used.\n" \
                    "Aborting 16S pipeline..."
        with pytest.raises(FileNotFoundError, match=re.escape(error_msg)):
            helpers.get_fastq_pass_dir(test_path)

    def test_fail_multiple_fastq_dirs(self):
        test_path = pathlib.Path(__file__).parent / "data" / "helpers" / "test_dir_multi_pass"
        error_msg = f"The directory {test_path} contains " \
                    f"multiple fastq_pass directories." \
                    "Please choose the one containing the fastq files you want to analyze and " \
                    "specify the entire path to the fastq_pass directory."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            helpers.get_fastq_pass_dir(test_path)
