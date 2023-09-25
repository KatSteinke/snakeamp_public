import pathlib
import re
import unittest

from unittest import mock

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


class TestFindExistingPath(unittest.TestCase):
    def test_find_all_exists(self):
        true_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_dir_2"\
                        / "rawdata"
        test_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_dir_2"\
                        / "rawdata"
        existing_path = snake_wrapper.get_existing_path(test_path)
        assert existing_path == true_path

    def test_find_partial(self):
        true_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_dir_4"
        test_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_dir_4" \
                    / "rawdata"
        existing_path = snake_wrapper.get_existing_path(test_path)
        assert existing_path == true_path


class TestSanitizePath(unittest.TestCase):
    def test_entire_path_exists_success(self):
        true_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_dir_2" \
                    / "rawdata"
        test_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_dir_2" \
                    / "rawdata"
        existing_path = snake_wrapper.get_clean_outdir(test_path)
        assert existing_path == true_path

    def test_space_in_existing(self):
        test_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test dir spaces"
        with pytest.raises(snake_wrapper.BadPathError,
                           match = "The path you are trying to save results to contains"
                                   " a space in an existing folder's name. "
                                   "This can break the pipeline. "
                                   "\nAborting...."):
            snake_wrapper.get_clean_outdir(test_path)

    @mock.patch(f'{snake_wrapper.__name__}.get_existing_path')
    def test_illegal_char_in_existing(self, mock_get_existing):
        mock_get_existing.return_value = "does_this_fail?"
        error_msg = "The path you are trying to save results to contains a character that " \
                    "can't be used in Windows in an existing folder's name. " \
                    "This can break the pipeline. " \
                    "\nAborting...."
        test_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_dir_2" \
                    / "rawdata"
        with pytest.raises(snake_wrapper.BadPathError, match = error_msg):
            snake_wrapper.get_clean_outdir(test_path)

    def test_reserved_name(self):
        plain_reserved = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_dir_2" \
                         / "rawdata" / "NUL"
        reserved_after_cleaning = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                                  / "test_dir_2" / "rawdata" / "N*UL"
        with pytest.raises(snake_wrapper.BadPathError,
                           match = "The path you are trying to save results to contains "
                                   "a name that is reserved in Windows. "
                                   "Cannot create this path. \n"
                                   "Aborting...."):
            snake_wrapper.get_clean_outdir(plain_reserved)
        with pytest.raises(snake_wrapper.BadPathError,
                           match = "The path you are trying to save results to contains "
                                   "a name that is reserved in Windows. "
                                   "Cannot create this path. \n"
                                   "Aborting...."):
            snake_wrapper.get_clean_outdir(reserved_after_cleaning)

    def test_strip_illegal_chars(self):
        messy_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_dir*5"
        cleaned_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_dir5"
        test_path = snake_wrapper.get_clean_outdir(messy_path)
        assert test_path == cleaned_path

    def test_strip_spaces(self):
        messy_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_dir 5"
        cleaned_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_dir_5"
        test_path = snake_wrapper.get_clean_outdir(messy_path)
        assert test_path == cleaned_path

    def test_new_path_success(self):
        clean_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_dir_5"
        test_path = snake_wrapper.get_clean_outdir(clean_path)
        assert test_path == clean_path


class TestGetExperimentName(unittest.TestCase):
    def test_experiment_name_success(self):
        runsheet = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                   / "runsheet_clean_name.xlsx"
        expected_name = "PLACEHOLDER_RUN_NAME"
        test_name = snake_wrapper.get_run_name(runsheet)
        assert expected_name == test_name

    def test_complain_blank_name(self):
        runsheet = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                   / "runsheet_blank_name.xlsx"
        error_msg = "No run name given in the runsheet."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            snake_wrapper.get_run_name(runsheet)

    def test_complain_no_name_column(self):
        runsheet = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                   / "runsheet_no_name_col.xlsx"
        error_msg = "No column giving the run name found in the runsheet."
        with pytest.raises(KeyError, match = re.escape(error_msg)):
            snake_wrapper.get_run_name(runsheet)


class TestGetNomadCommand(unittest.TestCase):
    def test_get_routine_command(self):
        indir = pathlib.Path("path/to/indir")
        outdir = pathlib.Path("path/to/outdir")
        runsheet =  pathlib.Path("path/to/runsheet")
        expected_command = ["nomad", "job", "dispatch",
                            "-meta", "indir=path/to/indir",
                            "-meta", "outdir=path/to/outdir",
                            "-meta", "runsheet=path/to/runsheet",
                            "16s-snake-emu-prod", snake_wrapper.default_config_file]
        test_command = snake_wrapper.get_pipeline_command(indir, outdir,runsheet)
        assert expected_command == test_command

    def test_get_test_command(self):
        """Test that the staging version is run if debug is specified (overriding default)."""
        indir = pathlib.Path("path/to/indir")
        outdir = pathlib.Path("path/to/outdir")
        runsheet = pathlib.Path("path/to/runsheet")
        expected_command = ["nomad", "job", "dispatch",
                            "-meta", "indir=path/to/indir",
                            "-meta", "outdir=path/to/outdir",
                            "-meta", "runsheet=path/to/runsheet",
                            "16s-snake-emu-staging", snake_wrapper.default_config_file]
        test_command = snake_wrapper.get_pipeline_command(indir, outdir, runsheet, debug = True)
        assert expected_command == test_command

    def test_run_different_config(self):
        """Test that a different configuration is used and given to the pipeline."""
        indir = pathlib.Path("path/to/indir")
        outdir = pathlib.Path("path/to/outdir")
        test_configfile = pathlib.Path("path/to/config")
        runsheet = pathlib.Path("path/to/runsheet")
        test_config = {"debug": True}
        expected_command = ["nomad", "job", "dispatch",
                            "-meta", "indir=path/to/indir",
                            "-meta", "outdir=path/to/outdir",
                            "-meta", "runsheet=path/to/runsheet",
                            "16s-snake-emu-staging", test_configfile]
        test_command = snake_wrapper.get_pipeline_command(indir, outdir, runsheet,
                                                          configfile = test_configfile,
                                                          active_config = test_config)
        assert expected_command == test_command

    def test_override_test_command(self):
        """Test that a config specifying debug mode can be overridden by the debug flag."""
        indir = pathlib.Path("path/to/indir")
        outdir = pathlib.Path("path/to/outdir")
        runsheet = pathlib.Path("path/to/runsheet")
        test_configfile = pathlib.Path("path/to/config")
        test_config = {"debug": True}
        expected_command = ["nomad", "job", "dispatch",
                            "-meta", "indir=path/to/indir",
                            "-meta", "outdir=path/to/outdir",
                            "-meta", "runsheet=path/to/runsheet",
                            "16s-snake-emu-prod", test_configfile]
        test_command = snake_wrapper.get_pipeline_command(indir, outdir, runsheet,
                                                          configfile = test_configfile,
                                                          active_config = test_config,
                                                          debug = False)
        assert expected_command == test_command
