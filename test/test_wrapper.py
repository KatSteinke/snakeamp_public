import pathlib
import re
import shutil
import unittest

from datetime import timedelta
from unittest import mock

import pandas as pd
import pytest

import monitor_run
import run_pipeline as snake_wrapper


class TestFindRundir(unittest.TestCase):
    def test_find_absolute_path_success(self):
        test_path = pathlib.Path(__file__).parent / "data" / "utilities_test" /"test_dir_2"
        minion_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "miniondir"
        true_path = test_path.resolve()
        pass_dirs = list(true_path.glob("rawdata/*/fastq_pass"))
        log_msg = "INFO:amplicon_nanopore:Data is retrieved from following folders: \n " \
                  f"{str([str(fastq_dir) for fastq_dir in pass_dirs])}"
        with self.assertLogs("amplicon_nanopore", level = "INFO") as logged:
            checked_path = snake_wrapper.find_rundir(test_path, minion_path)
            assert log_msg in logged.output
        assert checked_path == true_path

    def test_find_in_minion_dir_success(self):
        test_path = pathlib.Path("test1")
        minion_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "miniondir"
        true_path = pathlib.Path(__file__).parent / "data" / "utilities_test" /"miniondir" / "test1"
        pass_dirs = list(true_path.glob("rawdata/*/fastq_pass"))
        log_msg = "INFO:amplicon_nanopore:Data is retrieved from following folders: \n " \
                  f"{str([str(fastq_dir) for fastq_dir in pass_dirs])}"
        with self.assertLogs("amplicon_nanopore", level = "INFO") as logged:
            checked_path = snake_wrapper.find_rundir(test_path, minion_path)
            assert log_msg in logged.output
        assert checked_path == true_path

    def test_fail_path(self):
        test_path = pathlib.Path("test3")
        minion_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "miniondir"
        error_msg = "{} or {} does not exist \nAborting pipeline...".format(str(test_path),
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
Aborting pipeline...""".format(str(test_path))
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


class TestProcessRunsheet(unittest.TestCase):
    def test_catch_no_samples_for_analysis(self):
        """Fail if no samples have the amplicon type specified in the config."""
        runsheet = pathlib.Path(__file__).parent / "data" / "utilities_test"\
                   / "test_nanopore_runsheet.xlsx"
        active_config = {"sample_number_settings": {"sample_number_format":
                                                  '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})-\d?',
                                              "sample_numbers_in": "letter",
                                              "sample_numbers_out": "letter",
                                              "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                              "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                              "number_to_letter": {"70": "P",
                                                                   "30": "B",
                                                                   "10": "D",
                                                                   "50": "T"},
                                              "date_settings":
                                                  {"splice_in_date": False,
                                                   "length_without_date": 8,
                                                   "splice_after": 2},
                                              "negative_control": 'NegK[a-zA-Z0-9]*',
                                              "positive_control": {}},
                         "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                         "barcode_prefix": "RB",
                         "amplicon_type": "ITS"}
        error_msg = "No samples with amplicon type ITS found in runsheet"
        with pytest.raises(KeyError, match=re.escape(error_msg)):
            snake_wrapper.process_runsheet(runsheet, active_config)

    def test_filter_all_one_analysis(self):
        """Filter and check a runsheet in which all samples have the desired amplicon type."""
        runsheet = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                   / "test_nanopore_runsheet_16s_only.xlsx"
        active_config = {"sample_number_settings": {"sample_number_format":
                                                                      '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})-\d?',
                                                                  "sample_numbers_in": "letter",
                                                                  "sample_numbers_out": "letter",
                                                                  "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                                  "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                                  "number_to_letter": {"70": "P",
                                                                                       "30": "B",
                                                                                       "10": "D",
                                                                                       "50": "T"},
                                                                  "date_settings":
                                                                      {"splice_in_date": False,
                                                                       "length_without_date": 8,
                                                                       "splice_after": 2},
                                                                  "negative_control": 'NegK[a-zA-Z0-9]*',
                                                                  "positive_control": {}},
                                       "barcode_format": "NB[0-9]{2}",  # format of barcodes in runsheet
                                       "barcode_prefix": "NB",
                                       "amplicon_type": "16S"}
        expected_runsheet = pd.DataFrame(data={"Prøvenummer": ["F99123456-1",
                                                               "F99123456-2",
                                                               "NegK16S"],
                                               "Barkode": ["NB01", "NB42", "NB02"],
                                               "Eluat nr.": ["1", "2", "3"],
                                               "Analyse": ["16S", "16S", "16S"]})
        test_runsheet = snake_wrapper.process_runsheet(runsheet, active_config)
        pd.testing.assert_frame_equal(test_runsheet, expected_runsheet)

    def test_filter_multiple_types(self):
        """Filter and check a runsheet in which only some samples have the desired amplicon type."""
        runsheet = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                   / "test_nanopore_runsheet.xlsx"
        active_config = {"sample_number_settings": {"sample_number_format":
                                                        '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})-\d?',
                                                    "sample_numbers_in": "letter",
                                                    "sample_numbers_out": "letter",
                                                    "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                    "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                    "number_to_letter": {"70": "P",
                                                                         "30": "B",
                                                                         "10": "D",
                                                                         "50": "T"},
                                                    "date_settings":
                                                        {"splice_in_date": False,
                                                         "length_without_date": 8,
                                                         "splice_after": 2},
                                                    "negative_control": 'NegK[a-zA-Z0-9]*',
                                                    "positive_control": {}},
                         "barcode_format": "NB[0-9]{2}",  # format of barcodes in runsheet
                         "barcode_prefix": "NB",
                         "amplicon_type": "16S"}
        expected_runsheet = pd.DataFrame(data = {"Prøvenummer": ["F99123456-1",
                                                                 "F99123456-2",
                                                                 "NegK16S"],
                                                 "Barkode": ["NB01", "NB42", "NB02"],
                                                 "Eluat nr.": ["1", "2", "3"],
                                                 "Analyse": ["16S", "16S", "16S"]})
        test_runsheet = snake_wrapper.process_runsheet(runsheet, active_config)
        pd.testing.assert_frame_equal(test_runsheet, expected_runsheet)


class TestGetSeqTime(unittest.TestCase):
    default_time = 16

    @mock.patch("builtins.input")
    def test_fail_invalid_accept(self, mock_input):
        """Fail if an invalid value was entered for accepting/rejecting default sequencing time."""
        mock_input.return_value = "16"
        error_msg = "Sequencing time not entered. Aborting"
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            snake_wrapper.ask_seq_time(self.default_time)

    @mock.patch("builtins.input", side_effect=["n", "sixteen"])
    def test_fail_bad_time_format(self, mock_input):
        """Fail if a value that could not be converted to a time was entered."""
        error_msg = ("sixteen is not a valid sequencing time. "
                     "Sequencing time must be entered as numbers "
                     "(e.g. 8 for eight hours or 0.5 for half an hour).")
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            snake_wrapper.ask_seq_time(self.default_time)

    @mock.patch("builtins.input", side_effect = ["n", -1])
    def test_fail_invalid_time(self, mock_input):
        """Fail if an invalid time was entered."""
        error_msg = "Expected sequencing time must be greater than 0 hours."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            snake_wrapper.ask_seq_time(self.default_time)

    @mock.patch("builtins.input")
    def test_success_default_time(self, mock_input):
        """Return the default timespan if no changes are made."""
        mock_input.return_value = "y"
        expected_time = self.default_time
        test_time = snake_wrapper.ask_seq_time(self.default_time)
        assert expected_time == test_time

    @mock.patch("builtins.input", side_effect = ["n", 1])
    def test_success_different_time(self, mock_input):
        """Change the sequencing time when specified by the user."""
        expected_time = 1
        test_time = snake_wrapper.ask_seq_time(self.default_time)
        assert expected_time == test_time


class TestAskOutputPath(unittest.TestCase):
    default_path = pathlib.Path("data/test_run")

    @mock.patch("builtins.input")
    def test_fail_invalid_accept(self, mock_input):
        """Fail on an invalid response for whether to accept the default output path."""
        mock_input.return_value = "nope"
        error_msg = "Output folder not entered. Aborting"
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            snake_wrapper.ask_output_dir(self.default_path)

    # to check issues in the existing path, we need to pretend we've got a broken path
    # decorators are applied bottom up: https://stackoverflow.com/a/15922422/15704972
    @mock.patch(f'{snake_wrapper.__name__}.get_existing_path')
    @mock.patch("builtins.input", side_effect=["n", # user rejects path
                                               "/data/test:run"])  # user suggests new
    def test_fail_bad_user_path(self, mock_input, mock_path):
        """Fail when the user suggests a new path and it's invalid."""
        mock_path.return_value = "/data/test:run"
        error_msg = ("The path you are trying to save results to contains a character that "
                           "can't be used in Windows in an existing folder's name. "
                           "This can break the pipeline. "
                           "\nAborting....")
        with pytest.raises(snake_wrapper.BadPathError, match=re.escape(error_msg)):
            snake_wrapper.ask_output_dir(self.default_path)

    @mock.patch(f'{snake_wrapper.__name__}.get_existing_path')
    @mock.patch("builtins.input", side_effect=["/data/test:run/run2"])  # user needs to give a new path because the old one is broken
    def test_fail_bad_corrected_path(self, mock_input, mock_path):
        """Fail when the original path is invalid and the user's correction is as well."""
        mock_path.return_value = "/data/test:run"
        default_path = pathlib.Path("/data/test:run")
        error_msg = ("The path you are trying to save results to contains a character that "
                           "can't be used in Windows in an existing folder's name. "
                           "This can break the pipeline. "
                           "\nAborting....")
        log_msg = ("WARNING:amplicon_nanopore:The default target folder contains characters "
                   "that can break the pipeline.")
        with pytest.raises(snake_wrapper.BadPathError,
                           match=re.escape(error_msg)), self.assertLogs("amplicon_nanopore") as logged:
            snake_wrapper.ask_output_dir(default_path)
        assert log_msg in logged.output

    @mock.patch("builtins.input")
    def test_success_valid_default(self, mock_input):
        """Return the default output path when the user accepts it."""
        mock_input.return_value = "y"
        expected_path = self.default_path
        log_msg = f"INFO:amplicon_nanopore:Saving results to {self.default_path}"
        with self.assertLogs("amplicon_nanopore", level = "INFO") as logged:
            test_path = snake_wrapper.ask_output_dir(self.default_path)
            assert log_msg in logged.output
        assert expected_path == test_path

    @mock.patch("builtins.input", side_effect = ["n",  # user rejects the path
                                                 "/data/test_run2"])  # ...and gives a new one
    def test_success_valid_user_path(self, mock_input):
        """Return the user's new path if it's valid."""
        expected_path = pathlib.Path("/data/test_run2")
        log_msg = f"INFO:amplicon_nanopore:Saving results to {expected_path}"
        with self.assertLogs("amplicon_nanopore", level = "INFO") as logged:
            test_path = snake_wrapper.ask_output_dir(self.default_path)
            assert log_msg in logged.output
        assert expected_path == test_path

    @mock.patch("builtins.input", side_effect = ["/data/test_run2"])  # user needs to give a new path because the old one is broken
    def test_success_valid_correction(self, mock_input):
        """Ask the user for a new path and return it if the original path is invalid but the
        user's correction fixes it."""
        default_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test dir spaces"
        expected_path = pathlib.Path("/data/test_run2")
        log_msg = ("WARNING:amplicon_nanopore:The default target folder contains characters "
                   "that can break the pipeline.")
        success_msg = f"INFO:amplicon_nanopore:Saving results to {expected_path}"
        with self.assertLogs("amplicon_nanopore", level="INFO") as logged:
            test_path = snake_wrapper.ask_output_dir(default_path)
            assert log_msg in logged.output
            assert success_msg in logged.output
        assert expected_path == test_path

    @mock.patch("builtins.input")
    def test_success_correct_fixable_default(self, mock_input):
        """Correct a fixable bad path from default input."""
        mock_input.return_value = "y"
        default_path = pathlib.Path("data/test run")
        expected_path = self.default_path
        log_msg = f"INFO:amplicon_nanopore:Saving results to {expected_path}"
        with self.assertLogs("amplicon_nanopore", level = "INFO") as logged:
            test_path = snake_wrapper.ask_output_dir(default_path)
            assert log_msg in logged.output
        assert expected_path == test_path

    @mock.patch("builtins.input", side_effect = ["n",  # user rejects the path
                                                 "/data/test run2"])  # ...and gives a new one
    def test_success_correct_fixable_user(self, mock_input):
        """Correct a fixable bad path from user input."""
        expected_path = pathlib.Path("/data/test_run2")
        log_msg = f"INFO:amplicon_nanopore:Saving results to {expected_path}"
        with self.assertLogs("amplicon_nanopore") as logged:
            test_path = snake_wrapper.ask_output_dir(self.default_path)
            assert log_msg in logged.output
        assert expected_path == test_path


class TestCreateOutputDirs(unittest.TestCase):
    @classmethod
    def tearDownClass(cls) -> None:
        # clean up data after running
        shutil.rmtree(pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_output")
        shutil.rmtree(pathlib.Path(__file__).parent / "data" / "utilities_test"
                      / "test_existing_output" / "logs")

    def test_fail_existing_dir(self):
        """Fail if the output directory exists already."""
        error_msg = "The desired output directory already exists."
        with pytest.raises(FileExistsError, match = re.escape(error_msg)):
            snake_wrapper.set_up_output(pathlib.Path(__file__).parent / "data" / "utilities_test"
                                        / "test_existing_output")

    def test_success_create_new(self):
        """Create a new output directory and log directory."""
        output_dir = pathlib.Path(__file__).parent / "data" / "utilities_test"/ "test_output"
        assert not output_dir.exists()
        assert not (output_dir / "logs").exists()
        snake_wrapper.set_up_output(output_dir)
        assert output_dir.exists()
        assert (output_dir / "logs").exists()

    def test_success_continue_run(self):
        """Don't complain for a continued run, and log that it's continued."""
        log_msg = "INFO:amplicon_nanopore:Output directory already exists. Continuing run..."
        output_dir = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                      / "test_existing_output")
        assert output_dir.exists()
        assert not (output_dir / "logs").exists()
        with self.assertLogs("amplicon_nanopore", level = "INFO") as logged:
            snake_wrapper.set_up_output(output_dir, continue_run = True)
            assert log_msg in logged.output
        assert output_dir.exists()
        assert (output_dir / "logs").exists()


class TestInitializeRunFromInput(unittest.TestCase):
    active_config = {"sample_number_settings": {"sample_number_format":
                                                    '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})-\d?',
                                                "sample_numbers_in": "letter",
                                                "sample_numbers_out": "letter",
                                                "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                "number_to_letter": {"70": "P",
                                                                     "30": "B",
                                                                     "10": "D",
                                                                     "50": "T"},
                                                "date_settings":
                                                    {"splice_in_date": False,
                                                     "length_without_date": 8,
                                                     "splice_after": 2},
                                                "negative_control": 'NegK[a-zA-Z0-9]*',
                                                "positive_control": {}},
                     "barcode_format": "NB[0-9]{2}",  # format of barcodes in runsheet
                     "barcode_prefix": "NB",
                     "amplicon_type": "16S",
                     "debug": False,
                     "paths": {"output_base_path": "/path/to/output"},
                     "seq_run_duration_hours": 1}

    @mock.patch("builtins.input", side_effect = ["path/to/indir",  # sequencing directory
                                                 "y",  # accept sequencing time
                                                 str((pathlib.Path(__file__).parent / "data"  # runsheet
                                                      / "utilities_test"
                                                      / "test_nanopore_runsheet.xlsx")),
                                                 "y"])  # accept default outdir
    def test_success_use_default_time(self, mock_input):
        """Successfully set up a run using the default sequencing time."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        expected_outdir = pathlib.Path("/path/to/output/NANO_Amplicon_Y20990101_RUN0001_XYZ-16S")
        expected_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config,
                                               outdir = expected_outdir,
                                               sequencing_time = self.active_config["seq_run_duration_hours"])
        welcome_msg = ("INFO:amplicon_nanopore:### Nanopore 16S analysis\n"
                       "# Setup analysis -------------------------------")
        with self.assertLogs("amplicon_nanopore", level = "INFO") as logged:
            test_run = snake_wrapper.initialize_classic_run(active_config = self.active_config,
                                                            configfile = configfile)
            assert welcome_msg in logged.output
        assert test_run == expected_run

    @mock.patch("builtins.input", side_effect = ["path/to/indir",  # sequencing directory
                                                 "n",  # reject default sequencing time
                                                 1.5,  # set new sequencing time
                                                 str((pathlib.Path(
                                                     __file__).parent / "data"  # runsheet
                                                      / "utilities_test"
                                                      / "test_nanopore_runsheet.xlsx")),
                                                 "y"])  # accept default outdir
    def test_success_change_time(self, mock_input):
        """Successfully set up a run with a different sequencing time."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        expected_outdir = pathlib.Path("/path/to/output/NANO_Amplicon_Y20990101_RUN0001_XYZ-16S")
        expected_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config,
                                               outdir = expected_outdir,
                                               sequencing_time = 1.5)
        test_run = snake_wrapper.initialize_classic_run(active_config = self.active_config,
                                                        configfile = configfile)
        assert test_run == expected_run

    @mock.patch("builtins.input", side_effect = ["path/to/indir",  # sequencing directory
                                                "y",  # accept suggested sequencing time
                                                 str((pathlib.Path(
                                                     __file__).parent / "data"  # runsheet
                                                      / "utilities_test"
                                                      / "test_nanopore_runsheet.xlsx")),
                                                 "n",  # reject suggested outdir
                                                 "/path/to/new_outdir"])  # set new outdir
    def test_success_change_outdir(self, mock_input):
        """Successfully set up a run with a different output directory."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        expected_outdir = pathlib.Path("/path/to/new_outdir")
        expected_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config,
                                               outdir = expected_outdir,
                                               sequencing_time = self.active_config[
                                                   "seq_run_duration_hours"])
        test_run = snake_wrapper.initialize_classic_run(active_config = self.active_config,
                                                        configfile = configfile)
        assert test_run == expected_run

    @mock.patch("builtins.input", side_effect = ["path/to/indir",  # sequencing directory
                                                 "n",  # reject default sequencing time
                                                 1.5,  # set new sequencing time
                                                 str((pathlib.Path(
                                                     __file__).parent / "data"  # runsheet
                                                      / "utilities_test"
                                                      / "test_nanopore_runsheet.xlsx")),
                                                 "n",  # reject suggested outdir
                                                 "/path/to/new_outdir"])  # set new outdir
    def test_success_change_time_and_outdir(self, mock_input):
        """Successfully set up a run with a different sequencing time and output directory."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        expected_outdir = pathlib.Path("/path/to/new_outdir")
        expected_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config,
                                               outdir = expected_outdir,
                                               sequencing_time = 1.5)
        test_run = snake_wrapper.initialize_classic_run(active_config = self.active_config,
                                                        configfile = configfile)
        assert test_run == expected_run
