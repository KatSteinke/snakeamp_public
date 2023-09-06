import pathlib
import re
import unittest

import pytest

import run_pipeline as snake_wrapper


class TestFindRundir(unittest.TestCase):
    def test_find_absolute_path_success(self):
        test_path = pathlib.Path(__file__).parent / "data" / "utilities_test" /"test_dir_2"
        minion_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "miniondir"
        true_path = test_path.resolve()
        pass_dirs = list(true_path.glob("rawdata/*/fastq_pass"))
        log_msg = "INFO:16S_nanopore:Data is retrieved from following folders: \n " \
                  f"{str([str(fastq_dir) for fastq_dir in pass_dirs])}"
        with self.assertLogs("16S_nanopore", level = "INFO") as logged:
            checked_path = snake_wrapper.find_rundir(test_path, minion_path)
            assert log_msg in logged.output
        assert checked_path == true_path

    def test_find_in_minion_dir_success(self):
        test_path = pathlib.Path("test1")
        minion_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "miniondir"
        true_path = pathlib.Path(__file__).parent / "data" / "utilities_test" /"miniondir" / "test1"
        pass_dirs = list(true_path.glob("rawdata/*/fastq_pass"))
        log_msg = "INFO:16S_nanopore:Data is retrieved from following folders: \n " \
                  f"{str([str(fastq_dir) for fastq_dir in pass_dirs])}"
        with self.assertLogs("16S_nanopore", level = "INFO") as logged:
            checked_path = snake_wrapper.find_rundir(test_path, minion_path)
            assert log_msg in logged.output
        assert checked_path == true_path

    def test_fail_path(self):
        test_path = pathlib.Path("test3")
        minion_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "miniondir"
        error_msg = "{} or {} does not exist \nAborting ARTIC pipeline...".format(str(test_path),
                                                                                  str(minion_path
                                                                                      / "test3"))
        with pytest.raises(FileNotFoundError, match=re.escape(error_msg)):
            snake_wrapper.find_rundir(test_path, minion_path)

    def test_fail_rawdata(self):
        test_path = pathlib.Path(__file__).parent / "data" / "utilities_test" /"test_dir_4"
        minion_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "miniondir"
        error_msg = """fastq_pass folder(s) not found in expected location:
{}/rawdata/*/fastq_pass
Ensure correct directory and/or directory structure is used.
Aborting ARTIC pipeline...""".format(str(test_path))
        with pytest.raises(FileNotFoundError, match=re.escape(error_msg)):
            snake_wrapper.find_rundir(test_path, minion_path)
