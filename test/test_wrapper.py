import logging
import os
import pathlib
import re
import shutil
import unittest

from argparse import Namespace, ArgumentParser
from datetime import datetime, timedelta
from unittest import mock

import pandas as pd
import pytest

import monitor_run
import run_pipeline as snake_wrapper


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
                                                  r'([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})-\d?',
                                              "sample_numbers_in": "letter",
                                              "sample_numbers_out": "letter",
                                              "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                              "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                              "number_to_letter": {"70": "P",
                                                                   "30": "B",
                                                                   "10": "D",
                                                                   "11": "F",
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
                                                                      r'([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})-\d?',
                                                    "sample_numbers_in": "letter",
                                                    "sample_numbers_out": "letter",
                                                    "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                    "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                    "number_to_letter": {"70": "P",
                                                                         "30": "B",
                                                                         "10": "D",
                                                                         "11": "F",
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
                                                        r'([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})-\d?',
                                                    "sample_numbers_in": "letter",
                                                    "sample_numbers_out": "letter",
                                                    "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                    "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                    "number_to_letter": {"70": "P",
                                                                         "30": "B",
                                                                         "10": "D",
                                                                         "11": "F",
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

    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

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
        log_msg = ("The default target folder contains characters "
                   "that can break the pipeline.")
        with (pytest.raises(snake_wrapper.BadPathError,
                           match=re.escape(error_msg)),
              self._caplog.at_level(logging.WARNING, logger = "amplicon_nanopore")):
            snake_wrapper.ask_output_dir(default_path)
        assert ("amplicon_nanopore", logging.WARNING, log_msg) in self._caplog.record_tuples

    @mock.patch("builtins.input")
    def test_success_valid_default(self, mock_input):
        """Return the default output path when the user accepts it."""
        mock_input.return_value = "y"
        expected_path = self.default_path
        log_msg = f"Saving results to {self.default_path}"
        with self._caplog.at_level(logging.INFO, logger = "amplicon_nanopore"):
            test_path = snake_wrapper.ask_output_dir(self.default_path)
            assert ("amplicon_nanopore", logging.INFO, log_msg) in self._caplog.record_tuples
        assert expected_path == test_path

    @mock.patch("builtins.input", side_effect = ["n",  # user rejects the path
                                                 "/data/test_run2"])  # ...and gives a new one
    def test_success_valid_user_path(self, mock_input):
        """Return the user's new path if it's valid."""
        expected_path = pathlib.Path("/data/test_run2")
        log_msg = f"Saving results to {expected_path}"
        with self._caplog.at_level(logging.INFO, logger = "amplicon_nanopore"):
            test_path = snake_wrapper.ask_output_dir(self.default_path)
            assert ("amplicon_nanopore", logging.INFO, log_msg) in self._caplog.record_tuples
        assert expected_path == test_path

    @mock.patch("builtins.input", side_effect = ["/data/test_run2"])  # user needs to give a new path because the old one is broken
    def test_success_valid_correction(self, mock_input):
        """Ask the user for a new path and return it if the original path is invalid but the
        user's correction fixes it."""
        default_path = pathlib.Path(__file__).parent / "data" / "utilities_test" / "test dir spaces"
        expected_path = pathlib.Path("/data/test_run2")
        log_msg = ("The default target folder contains characters "
                   "that can break the pipeline.")
        success_msg = f"Saving results to {expected_path}"
        with self._caplog.at_level(logging.INFO, logger = "amplicon_nanopore"):
            test_path = snake_wrapper.ask_output_dir(default_path)
            assert ("amplicon_nanopore", logging.WARNING, log_msg) in self._caplog.record_tuples
            assert ("amplicon_nanopore", logging.INFO, success_msg) in self._caplog.record_tuples
        assert expected_path == test_path

    @mock.patch("builtins.input")
    def test_success_correct_fixable_default(self, mock_input):
        """Correct a fixable bad path from default input."""
        mock_input.return_value = "y"
        default_path = pathlib.Path("data/test run")
        expected_path = self.default_path
        log_msg = f"Saving results to {expected_path}"
        with self._caplog.at_level(logging.INFO, logger = "amplicon_nanopore"):
            test_path = snake_wrapper.ask_output_dir(default_path)
            assert ("amplicon_nanopore", logging.INFO, log_msg) in self._caplog.record_tuples
        assert expected_path == test_path

    @mock.patch("builtins.input", side_effect = ["n",  # user rejects the path
                                                 "/data/test run2"])  # ...and gives a new one
    def test_success_correct_fixable_user(self, mock_input):
        """Correct a fixable bad path from user input."""
        expected_path = pathlib.Path("/data/test_run2")
        log_msg = f"Saving results to {expected_path}"  # TODO: be more explicit about the path being changed?
        with self._caplog.at_level(logging.INFO, logger = "amplicon_nanopore"):
            test_path = snake_wrapper.ask_output_dir(self.default_path)
            assert ("amplicon_nanopore", logging.INFO, log_msg) in self._caplog.record_tuples
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
                                                    r'([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})-\d?',
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

    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    @mock.patch("builtins.input")
    def test_success_use_default_time(self, mock_input):
        """Successfully set up a run using the default sequencing time."""
        expected_indir = (pathlib.Path(__file__).parent / "data"/"utilities_test"/"miniondir"
                          /"test1"/"rawdata"/"test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        expected_outdir = pathlib.Path("/path/to/output/NANO_Amplicon_Y20990101_RUN0001_XYZ-16S")
        mock_input.side_effect = [str(expected_indir),  # sequencing directory
                                  "y",  # accept sequencing time
                                  str(runsheet),  # runsheet
                                  "y"]  # accept default outdir
        expected_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config,
                                               outdir = expected_outdir,
                                               sequencing_time = self.active_config[
                                                   "seq_run_duration_hours"])
        welcome_msg = ("### Nanopore 16S analysis\n"
                       "# Setup analysis -------------------------------")
        with self._caplog.at_level(logging.INFO, logger = "amplicon_nanopore"):
            test_run = snake_wrapper.initialize_classic_run(active_config = self.active_config,
                                                            configfile = configfile)
            assert ("amplicon_nanopore", logging.INFO, welcome_msg) in self._caplog.record_tuples
        assert test_run == expected_run

    @mock.patch("builtins.input")
    def test_success_test_run(self, mock_input):
        """Successfully set up a run set as a test run in the config."""
        expected_indir = (pathlib.Path(__file__).parent / "data"/"utilities_test"/"miniondir"
                          /"test1"/"rawdata"/"test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        # set debug param
        active_config = self.active_config.copy()
        active_config["debug"] = True
        expected_outdir = pathlib.Path("/path/to/output/NANO_Amplicon_Y20990101_RUN0001_XYZ-16S")
        mock_input.side_effect = [str(expected_indir),  # sequencing directory
                                  "y",  # accept sequencing time
                                  str(runsheet),  # runsheet
                                  "y"]  # accept default outdir
        expected_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = active_config,
                                               outdir = expected_outdir,
                                               sequencing_time = active_config[
                                                   "seq_run_duration_hours"],
                                               test_run = True)
        welcome_msg = ("### Nanopore 16S analysis\n"
                       "# Setup analysis -------------------------------")
        with self._caplog.at_level(logging.INFO, logger = "amplicon_nanopore"):
            test_run = snake_wrapper.initialize_classic_run(active_config = active_config,
                                                            configfile = configfile)
            assert ("amplicon_nanopore", logging.INFO, welcome_msg) in self._caplog.record_tuples
        assert test_run == expected_run

    @mock.patch("builtins.input")
    def test_success_change_time(self, mock_input):
        """Successfully set up a run with a different sequencing time."""
        expected_indir = (pathlib.Path(__file__).parent / "data"/"utilities_test"/"miniondir"
                          /"test1"/"rawdata"/"test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        expected_outdir = pathlib.Path("/path/to/output/NANO_Amplicon_Y20990101_RUN0001_XYZ-16S")
        mock_input.side_effect = [str(expected_indir),  # sequencing directory
                                  "n",  # reject default sequencing time
                                  1.5,  # set new sequencing time
                                  str(runsheet),  # runsheet
                                  "y"]  # accept default outdir
        expected_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config,
                                               outdir = expected_outdir, sequencing_time = 1.5)
        test_run = snake_wrapper.initialize_classic_run(active_config = self.active_config,
                                                        configfile = configfile)
        assert test_run == expected_run

    @mock.patch("builtins.input")
    def test_success_change_outdir(self, mock_input):
        """Successfully set up a run with a different output directory."""
        expected_indir = (pathlib.Path(__file__).parent / "data"/"utilities_test"/"miniondir"
                          /"test1"/"rawdata"/"test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        expected_outdir = pathlib.Path("/path/to/new_outdir")
        mock_input.side_effect = [str(expected_indir),  # sequencing directory
                                  "y",  # accept suggested sequencing time
                                  str(runsheet),
                                  "n",  # reject suggested outdir
                                  "/path/to/new_outdir"]  # set new outdir
        expected_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config,
                                               outdir = expected_outdir,
                                               sequencing_time = self.active_config[
                                                   "seq_run_duration_hours"])
        test_run = snake_wrapper.initialize_classic_run(active_config = self.active_config,
                                                        configfile = configfile)
        assert test_run == expected_run

    @mock.patch("builtins.input")
    def test_success_change_time_and_outdir(self, mock_input):
        """Successfully set up a run with a different sequencing time and output directory."""
        expected_indir = (pathlib.Path(__file__).parent / "data"/"utilities_test"/"miniondir"
                          /"test1"/"rawdata"/"test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        expected_outdir = pathlib.Path("/path/to/new_outdir")
        mock_input.side_effect = [str(expected_indir),  # sequencing directory
                                  "n",  # reject default sequencing time
                                  1.5,  # set new sequencing time
                                  str(runsheet),  # runsheet
                                  "n",  # reject suggested outdir
                                  "/path/to/new_outdir"]  # set new outdir
        expected_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config,
                                               outdir = expected_outdir, sequencing_time = 1.5)
        test_run = snake_wrapper.initialize_classic_run(active_config = self.active_config,
                                                        configfile = configfile)
        assert test_run == expected_run
class TestInitializeRunFromCommandline(unittest.TestCase):
    active_config = {"input_names": pathlib.Path(__file__).parent / "data" / "input_names"
                                    / "input_da_old_lis.yaml",
                     "sample_number_settings": {"sample_number_format":
                                                    r'([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})-\d?',
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
                     "seq_run_duration_hours": 1,
                     "lab_info_system":
                         {"use_lis_features": True,
                          "lis_report": str(pathlib.Path(__file__).parent / "data" / "summarize_emu"
                                            / "fake_mads_material.csv"),
                          "database": "",
                          "dialect": ""}
                     }
    fetch_lis_config = {"input_names": pathlib.Path(__file__).parent / "data" / "input_names"
                                       / "input_da_old_lis.yaml",
                        "sample_number_settings": {"sample_number_format":
                                                       r'([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})-\d?',
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
                        "seq_run_duration_hours": 1,
                        "lab_info_system":
                            {"use_lis_features": True,
                             "lis_report": str(pathlib.Path(__file__).parent
                                               / "data"
                                               / "utilities_test"
                                               / "test_lis_from_db.txt"),
                             "database": "/:memory:",
                             "dialect": "sqlite+pysqlite",
                             "filter_column": "isol_16s",
                             "filter_value": "16S",
                             "db_table": "test_db",
                             "schema": ""}
                        }
    indir = (pathlib.Path(__file__).parent / "data" / "utilities_test" / "miniondir"
             / "test1" / "rawdata" / "test_subdir")
    runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                / "test_nanopore_runsheet.xlsx")
    configfile = pathlib.Path("path/to/config")
    expected_outdir = pathlib.Path("/path/to/output"
                                   "/NANO_Amplicon_Y20990101_RUN0001_XYZ-16S")

    def tearDown(self) -> None:
        """Make sure to blank out the text of the output file."""
        (pathlib.Path(__file__).parent
            / "data"
            / "utilities_test"
            / "test_lis_from_db.txt").write_text(data = "", encoding = "latin-1")

    def test_missing_runsheet(self):
        """Fail if runsheet was not entered."""
        args = Namespace(runsheet = None, rundir = self.indir, test_run=None, run_time = None,
                         outdir = None)
        error_msg = "Runsheet not specified."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            snake_wrapper.initialize_commandline_run(args)

    def test_missing_data(self):
        """Fail if path to sequencing data was not entered."""
        args = Namespace(runsheet = self.runsheet, rundir = None, test_run = False, run_time = None,
                         outdir = None)
        error_msg = "Nanopore run directory not specified."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            snake_wrapper.initialize_commandline_run(args)

    # overwrite the actual default config so we have something stable
    @mock.patch.dict(f"{snake_wrapper.__name__}.WORKFLOW_CONFIG", active_config,
                     clear = True)
    def test_default_run(self):
        """Set up a run with default settings."""
        args = Namespace(runsheet = self.runsheet, rundir = self.indir, test_run = None,
                         run_time = None, outdir = None)
        expected_run = monitor_run.AmpliconRun(sequence_dir = self.indir, runsheet = self.runsheet,
                                               configfile = snake_wrapper.DEFAULT_CONFIG_FILE,
                                               active_config = self.active_config,
                                               outdir = self.expected_outdir,
                                               sequencing_time = self.active_config[
                                                   "seq_run_duration_hours"])
        test_run = snake_wrapper.initialize_commandline_run(args)
        assert expected_run == test_run

    def test_set_config(self):
        """Set up a run with settings set through a different config."""
        args = Namespace(runsheet = self.runsheet, rundir = self.indir, test_run = None,
                         run_time = None, outdir = None)
        test_config = self.active_config.copy()
        test_config["seq_run_duration_hours"] = 0.5
        expected_run = monitor_run.AmpliconRun(sequence_dir = self.indir, runsheet = self.runsheet,
                                               configfile = self.configfile,
                                               active_config = test_config,
                                               outdir = self.expected_outdir, sequencing_time = 0.5)
        test_run = snake_wrapper.initialize_commandline_run(args, test_config, self.configfile)
        assert expected_run == test_run

    def test_no_lis(self):
        """Run correctly even if no LIS is specified."""
        args = Namespace(runsheet = self.runsheet, rundir = self.indir, test_run = None,
                         run_time = None, outdir = None)
        test_config = self.active_config.copy()
        test_config["lab_info_system"]["use_lis"] = False
        test_config["lab_info_system"]["lis_report"] = ""
        expected_run = monitor_run.AmpliconRun(sequence_dir = self.indir, runsheet = self.runsheet,
                                               configfile = self.configfile,
                                               active_config = test_config,
                                               outdir = self.expected_outdir, sequencing_time = 1)
        test_run = snake_wrapper.initialize_commandline_run(args, test_config, self.configfile)
        assert expected_run == test_run

    def test_set_run_time(self):
        """Override run time from the commandline."""
        args = Namespace(runsheet = self.runsheet, rundir = self.indir, test_run = None,
                         run_time = 0.5, outdir = None)
        expected_run = monitor_run.AmpliconRun(sequence_dir = self.indir, runsheet = self.runsheet,
                                               configfile = self.configfile,
                                               active_config = self.active_config,
                                               outdir = self.expected_outdir, sequencing_time = 0.5)
        test_run = snake_wrapper.initialize_commandline_run(args, self.active_config,
                                                            self.configfile)
        assert expected_run == test_run

    def test_set_outdir(self):
        """Override output directory from the commandline."""
        args = Namespace(runsheet = self.runsheet, rundir = self.indir, test_run = None,
                         run_time = None, outdir = "path/to/test_outdir")
        expected_run = monitor_run.AmpliconRun(sequence_dir = self.indir, runsheet = self.runsheet,
                                               configfile = self.configfile,
                                               active_config = self.active_config,
                                               outdir = pathlib.Path("path/to/test_outdir"),
                                               sequencing_time = self.active_config[
                                                   "seq_run_duration_hours"])
        test_run = snake_wrapper.initialize_commandline_run(args, self.active_config,
                                                            self.configfile)
        assert expected_run == test_run

    def test_set_test_mode(self):
        """Set debug mode from the commandline."""
        args = Namespace(runsheet = self.runsheet, rundir = self.indir, test_run = True,
                         run_time = None, outdir = None)
        expected_run = monitor_run.AmpliconRun(sequence_dir = self.indir, runsheet = self.runsheet,
                                               configfile = self.configfile,
                                               active_config = self.active_config,
                                               outdir = self.expected_outdir,
                                               sequencing_time = self.active_config[
                                                   "seq_run_duration_hours"], test_run = True)
        test_run = snake_wrapper.initialize_commandline_run(args, self.active_config,
                                                            self.configfile)
        assert expected_run == test_run


class TestCheckClassicMode(unittest.TestCase):
    def test_fail_extra_start_args(self):
        """Fail if the pipeline is started with arguments that aren't in the parent parser."""
        fail_args = Namespace(foo = None, bar = None)
        parser = ArgumentParser(add_help = False)
        parser.add_argument("--baz")
        error_msg = "Argument(s) ['bar', 'foo'] are not valid arguments for the specified parser."
        with pytest.raises(KeyError, match = re.escape(error_msg)):
            snake_wrapper.check_if_classic_mode(fail_args, parser)

    def test_fail_extra_classic_mode_args(self):
        """Fail if the permitted classic mode args contain args not in the parent parser."""
        true_args = Namespace(baz=None)
        fail_args = ["foo", "bar"]
        parser = ArgumentParser(add_help = False)
        parser.add_argument("--baz")
        error_msg = "Argument(s) ['bar', 'foo'] are not valid arguments for the specified parser."
        with pytest.raises(KeyError, match = re.escape(error_msg)):
            snake_wrapper.check_if_classic_mode(true_args, parser, classic_mode_args = fail_args)

    def test_success_classic_mode_all_falsy(self):
        """Successfully recognize classic mode when everything is None/False."""
        test_args = Namespace(foo = None, bar = None)
        parser = ArgumentParser(add_help = False)
        parser.add_argument("--foo")
        parser.add_argument("--bar")
        assert snake_wrapper.check_if_classic_mode(test_args, parser) is True

    def test_success_classic_mode_some_truthy(self):
        """Successfully recognize classic mode when some default arguments default to truthy
         values."""
        test_args = Namespace(foo = 42, bar = False)
        parser = ArgumentParser(add_help = False)
        parser.add_argument("--foo", default = 42)
        parser.add_argument("--bar", action = "store_true")
        assert snake_wrapper.check_if_classic_mode(test_args, parser) is True

    def test_success_store_false(self):
        """Report parser called without args when the parser has an argument with store_false
        ( = called without the flag = True)."""
        test_args = Namespace(foo = None, bar = True)
        parser = ArgumentParser(add_help = False)
        parser.add_argument("--foo")
        parser.add_argument("--bar", action = "store_false")
        assert snake_wrapper.check_if_classic_mode(test_args, parser) is True

    def test_success_classic_mode_args_set(self):
        """Successfully recognize classic mode when some parameters can be set."""
        test_args = Namespace(foo = 42, bar = False)
        classic_mode_args = ["foo"]
        parser = ArgumentParser(add_help = False)
        parser.add_argument("--foo")
        parser.add_argument("--bar", action = "store_true")
        assert snake_wrapper.check_if_classic_mode(test_args, parser, classic_mode_args) is True

    def test_success_not_classic_mode(self):
        """Successfully recognize non-default mode."""
        test_args = Namespace(foo = 42, bar = None)
        parser = ArgumentParser(add_help = False)
        parser.add_argument("--foo")
        parser.add_argument("--bar")
        assert snake_wrapper.check_if_classic_mode(test_args, parser) is False


class TestRunPipeline(unittest.TestCase):
    active_config = {"input_names": pathlib.Path(__file__).parent.resolve() / "data" / "input_names"
                                    / "input_da_old_lis.yaml",
                     "sample_number_settings": {"sample_number_format":
                                                    r'([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})-\d?',
                                                "sample_numbers_in": "letter",
                                                "sample_numbers_out": "letter",
                                                "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                "number_to_letter": {"70": "P",
                                                                     "30": "B",
                                                                     "10": "D",
                                                                     "11": "F",
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
                     "debug": True,
                     "run_on": "nomad",
                     "cores": 8,
                     "paths": {"output_base_path": str(pathlib.Path(__file__).parent / "data"
                                                       / "utilities_test" / "test_outdir")},
                     "seq_run_duration_hours": 1,
                     "check_interval_seconds": 1,
                     "lab_info_system":
                         {"use_lis_features": True,
                          "lis_report": str(pathlib.Path(__file__).parent / "data" / "summarize_emu"
                                            / "fake_mads_material.csv"),
                          "database": "",
                          "dialect": ""}
                     }
    active_config_file = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                          / "test_16s_config.yaml").resolve()

    def tearDown(self):
        # clean up any existing paths
        logfiles = [pathlib.Path(__file__).parent / "data" / "utilities_test" / "new_logfile.log",
                    (pathlib.Path(__file__).parent / "data" / "utilities_test"
                     / "test_existing_output" / "logs" / "start_pipeline.log")
                    ]
        for logfile in logfiles:
            if logfile.exists():
                os.unlink(logfile)

        output_paths = [(pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_outdir"
                         / "NANO_Amplicon_Y20990101_RUN0001_XYZ-16S"),
                        (pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_outdir"
                         / "NANO_Amplicon_Y20990101_RUN0001_XYZ-18S"),
                        (pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_outdir"
                         / "test_manual_outdir_classic"),
                        (pathlib.Path(__file__).parent / "data" / "utilities_test" / "test_outdir"
                         / "test_manual_outdir_commandline"),
                        (pathlib.Path(__file__).parent / "data" / "utilities_test"
                         / "test_existing_output" / "logs")
                        ]
        for path in output_paths:
            if path.exists():
                shutil.rmtree(path)


        mads_db_file = (pathlib.Path(__file__).parent
                        / "data"
                        / "utilities_test"
                        / "test_lis_from_db.txt")
        with open(mads_db_file, "w") as mads_db:
            mads_db.write("")

    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    # mock default config: something reusable
    # TODO shouldn't this fail if the check interval is missing...
    @mock.patch.dict(f"{snake_wrapper.__name__}.WORKFLOW_CONFIG", active_config,
                     clear = True)
    @mock.patch(f"{snake_wrapper.__name__}.DEFAULT_CONFIG_FILE", active_config_file)
    # mock classic mode check so we can get away with using dry run
    @mock.patch(f"{snake_wrapper.__name__}.check_if_classic_mode", return_value=True)
    # mock fork so it doesn't actually fork off anything
    @mock.patch(f"{snake_wrapper.__name__}.os.fork")
    # mock input
    @mock.patch("builtins.input")
    # TODO: set up logging to file and test that
    def test_run_classic_default(self, mock_input, mock_fork, mock_mode):
        """Start a run in classic mode with the default config."""
        mock_fork.return_value = False
        expected_indir = (pathlib.Path(__file__).parent / "data" / "monitor_run" / "miniondir"
                          / "test1" / "rawdata" / "test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = self.active_config_file
        expected_outdir = pathlib.Path(self.active_config["paths"]["output_base_path"]) \
                                       / "NANO_Amplicon_Y20990101_RUN0001_XYZ-16S"
        mock_input.side_effect = [str(expected_indir),  # sequencing directory
                                  "y",  # accept sequencing time
                                  str(runsheet),  # runsheet
                                  "y"]  # accept default outdir
        assert not expected_outdir.exists()
        nomad_command = ["nomad", "job", "dispatch",
                         "-meta", f"indir={expected_indir}",
                         "-meta", f"outdir={expected_outdir}",
                         "-meta", f"runsheet={runsheet}",
                         "16s-snake-emu-staging", str(configfile)]
        expected_command = ["echo", f'"{" ".join(nomad_command)}"']
        test_command = snake_wrapper.run_pipeline(["--dry_run"])
        assert expected_command == test_command.args
        assert expected_outdir.exists()

    # TODO: test logging
    @mock.patch(f"{snake_wrapper.__name__}.check_if_classic_mode", return_value=True)
    # mock fork so it doesn't actually fork off anything
    @mock.patch(f"{snake_wrapper.__name__}.os.fork")
    # mock input
    @mock.patch("builtins.input")
    def test_run_classic_set_config(self, mock_input, mock_fork, mock_mode):
        """Start a run in classic mode while giving a different config."""
        mock_fork.return_value = False
        expected_indir = (pathlib.Path(__file__).parent / "data" / "monitor_run" / "miniondir"
                          / "test1" / "rawdata" / "test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = str(pathlib.Path(__file__).parent / "data"
                                                   /"utilities_test"/"test_18s_config.yaml")
        expected_outdir = (pathlib.Path(self.active_config["paths"]["output_base_path"]).
                           relative_to(pathlib.Path(__file__).parent.parent)
                          / "NANO_Amplicon_Y20990101_RUN0001_XYZ-18S")
        mock_input.side_effect = [str(expected_indir),  # sequencing directory
                                  "y",  # accept sequencing time
                                  str(runsheet),  # runsheet
                                  "y"]  # accept default outdir
        assert not expected_outdir.exists()
        nomad_command = ["nomad", "job", "dispatch",
                         "-meta", f"indir={expected_indir}",
                         "-meta", f"outdir={expected_outdir}",
                         "-meta", f"runsheet={runsheet}",
                         "16s-snake-emu-staging", str(configfile)]
        expected_command = ["echo", f'"{" ".join(nomad_command)}"']
        test_args = ["--workflow_config_file", str(configfile), "--dry_run"]
        test_command = snake_wrapper.run_pipeline(test_args)
        assert expected_command == test_command.args
        assert expected_outdir.exists()

    @mock.patch(f"{snake_wrapper.__name__}.check_if_classic_mode", return_value=True)
    # mock fork so it doesn't actually fork off anything
    @mock.patch(f"{snake_wrapper.__name__}.os.fork")
    # mock input
    @mock.patch("builtins.input")
    def test_run_classic_local(self, mock_input, mock_fork, mock_mode):
        """Start a local run in classic mode."""
        mock_fork.return_value = False
        expected_indir = (pathlib.Path(__file__).parent / "data" / "monitor_run" / "miniondir"
                          / "test1" / "rawdata" / "test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = str(pathlib.Path(__file__).parent / "data"
                                                   /"utilities_test"/"test_18s_config_local.yaml")
        expected_outdir = (pathlib.Path(self.active_config["paths"]["output_base_path"]).
                           relative_to(pathlib.Path(__file__).parent.parent)
                          / "NANO_Amplicon_Y20990101_RUN0001_XYZ-18S")
        mock_input.side_effect = [str(expected_indir),  # sequencing directory
                                  "y",  # accept sequencing time
                                  str(runsheet),  # runsheet
                                  "y"]  # accept default outdir
        assert not expected_outdir.exists()
        local_command = ["snakemake", "-s", "Snakefile",
                            "--cores", "8",
                            "--keep-going",
                            "--config",
                            f"outdir={expected_outdir}",
                            f"rundir={expected_indir}",
                            f"runsheet={runsheet}",
                            f"config_path={configfile}",
                            "--configfile", configfile]
        expected_command = ["echo", f'"{" ".join(local_command)}"']
        test_args = ["--workflow_config_file", str(configfile), "--dry_run"]
        test_command = snake_wrapper.run_pipeline(test_args)
        assert expected_command == test_command.args
        assert expected_outdir.exists()

    @mock.patch.dict(f"{snake_wrapper.__name__}.WORKFLOW_CONFIG", active_config,
                     clear = True)
    @mock.patch(f"{snake_wrapper.__name__}.DEFAULT_CONFIG_FILE", active_config_file)
    # mock classic mode check so we can get away with using dry run
    @mock.patch(f"{snake_wrapper.__name__}.check_if_classic_mode", return_value = True)
    # mock fork so it doesn't actually fork off anything
    @mock.patch(f"{snake_wrapper.__name__}.os.fork")
    # mock input
    @mock.patch("builtins.input")
    def test_set_manual_outdir_classic(self, mock_input, mock_fork, mock_mode):
        """Start a run in classic mode and set the output directory."""
        mock_fork.return_value = False
        expected_indir = (pathlib.Path(__file__).parent / "data" / "monitor_run" / "miniondir"
                          / "test1" / "rawdata" / "test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = self.active_config_file
        expected_outdir = pathlib.Path(self.active_config["paths"]["output_base_path"]) \
                          / "NANO_Amplicon_Y20990101_RUN0001_XYZ-16S"
        mock_input.side_effect = [str(expected_indir),  # sequencing directory
                                  "y",  # accept sequencing time
                                  str(runsheet),  # runsheet
                                  "n", # don't accept default outdir
                                  str(expected_outdir)]
        assert not expected_outdir.exists()
        nomad_command = ["nomad", "job", "dispatch",
                         "-meta", f"indir={expected_indir}",
                         "-meta", f"outdir={expected_outdir}",
                         "-meta", f"runsheet={runsheet}",
                         "16s-snake-emu-staging", str(configfile)]
        expected_command = ["echo", f'"{" ".join(nomad_command)}"']
        test_command = snake_wrapper.run_pipeline(["--dry_run"])
        assert expected_command == test_command.args
        assert expected_outdir.exists()

    @mock.patch.dict(f"{snake_wrapper.__name__}.WORKFLOW_CONFIG", active_config,
                     clear = True)
    @mock.patch(f"{snake_wrapper.__name__}.DEFAULT_CONFIG_FILE", active_config_file)
    # mock fork so it doesn't actually fork off anything
    @mock.patch(f"{snake_wrapper.__name__}.os.fork")
    def test_run_commandline_default(self, mock_fork):
        """Start a run in commandline mode with the default config."""
        mock_fork.return_value = False
        expected_indir = (pathlib.Path(__file__).parent / "data" / "monitor_run" / "miniondir"
                          / "test1" / "rawdata" / "test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = self.active_config_file
        expected_outdir = pathlib.Path(self.active_config["paths"]["output_base_path"]) \
                          / "NANO_Amplicon_Y20990101_RUN0001_XYZ-16S"
        assert not expected_outdir.exists()
        nomad_command = ["nomad", "job", "dispatch",
                         "-meta", f"indir={expected_indir}",
                         "-meta", f"outdir={expected_outdir}",
                         "-meta", f"runsheet={runsheet}",
                         "16s-snake-emu-staging", str(configfile)]
        expected_command = ["echo", f'"{" ".join(nomad_command)}"']
        test_args = ["--rundir", str(expected_indir), "--runsheet", str(runsheet), "--dry_run"]
        test_command = snake_wrapper.run_pipeline(test_args)
        assert expected_command == test_command.args
        assert expected_outdir.exists()

    @mock.patch(f"{snake_wrapper.__name__}.os.fork")
    def test_run_commandline_set_config(self, mock_fork):
        """Start a run in commandline mode while giving a different config."""
        mock_fork.return_value = False
        expected_indir = (pathlib.Path(__file__).parent / "data" / "monitor_run" / "miniondir"
                          / "test1" / "rawdata" / "test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = str(pathlib.Path(__file__).parent / "data"
                         / "utilities_test" / "test_18s_config.yaml")
        expected_outdir = (pathlib.Path(self.active_config["paths"]["output_base_path"]).
                           relative_to(pathlib.Path(__file__).parent.parent)
                          / "NANO_Amplicon_Y20990101_RUN0001_XYZ-18S")
        assert not expected_outdir.exists()
        nomad_command = ["nomad", "job", "dispatch",
                         "-meta", f"indir={expected_indir}",
                         "-meta", f"outdir={expected_outdir}",
                         "-meta", f"runsheet={runsheet}",
                         "16s-snake-emu-staging", str(configfile)]
        expected_command = ["echo", f'"{" ".join(nomad_command)}"']
        test_args = ["--rundir", str(expected_indir), "--runsheet", str(runsheet),
                     "--workflow_config_file", str(configfile), "--dry_run"]
        test_command = snake_wrapper.run_pipeline(test_args)
        assert expected_command == test_command.args
        assert expected_outdir.exists()

    @mock.patch(f"{snake_wrapper.__name__}.os.fork")
    def test_run_commandline_local(self, mock_fork):
        """Start a local run in commandline mode."""
        mock_fork.return_value = False
        expected_indir = (pathlib.Path(__file__).parent / "data" / "monitor_run" / "miniondir"
                          / "test1" / "rawdata" / "test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = str(pathlib.Path(__file__).parent / "data"
                         / "utilities_test" / "test_18s_config_local.yaml")
        expected_outdir = (pathlib.Path(self.active_config["paths"]["output_base_path"]).
                           relative_to(pathlib.Path(__file__).parent.parent)
                          / "NANO_Amplicon_Y20990101_RUN0001_XYZ-18S")
        assert not expected_outdir.exists()
        local_command = ["snakemake", "-s", "Snakefile",
                            "--cores", "8",
                            "--keep-going",
                            "--config",
                            f"outdir={expected_outdir}",
                            f"rundir={expected_indir}",
                            f"runsheet={runsheet}",
                            f"config_path={configfile}",
                            "--configfile", configfile]
        expected_command = ["echo", f'"{" ".join(local_command)}"']
        test_args = ["--rundir", str(expected_indir), "--runsheet", str(runsheet),
                     "--workflow_config_file", str(configfile), "--dry_run"]
        test_command = snake_wrapper.run_pipeline(test_args)
        assert expected_command == test_command.args
        assert expected_outdir.exists()

    @mock.patch.dict(f"{snake_wrapper.__name__}.WORKFLOW_CONFIG", active_config,
                     clear = True)
    @mock.patch(f"{snake_wrapper.__name__}.DEFAULT_CONFIG_FILE", active_config_file)
    # mock fork so it doesn't actually fork off anything
    @mock.patch(f"{snake_wrapper.__name__}.os.fork")
    def test_run_commandline_set_logfile(self, mock_fork):
        """Start a run in commandline mode and specify a logfile."""
        mock_fork.return_value = False
        expected_indir = (pathlib.Path(__file__).parent / "data" / "monitor_run" / "miniondir"
                          / "test1" / "rawdata" / "test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = self.active_config_file
        expected_outdir = pathlib.Path(self.active_config["paths"]["output_base_path"]) \
                          / "NANO_Amplicon_Y20990101_RUN0001_XYZ-16S"
        new_logfile = pathlib.Path(__file__).parent / "data" / "utilities_test" / "new_logfile.log"
        assert not expected_outdir.exists()
        assert not new_logfile.exists()
        nomad_command = ["nomad", "job", "dispatch",
                         "-meta", f"indir={expected_indir}",
                         "-meta", f"outdir={expected_outdir}",
                         "-meta", f"runsheet={runsheet}",
                         "16s-snake-emu-staging", str(configfile)]
        expected_command = ["echo", f'"{" ".join(nomad_command)}"']
        test_args = ["--rundir", str(expected_indir), "--runsheet", str(runsheet),
                     "--logfile", str(new_logfile),
                     "--dry_run"]
        test_command = snake_wrapper.run_pipeline(test_args)
        assert expected_command == test_command.args
        assert expected_outdir.exists()
        assert new_logfile.exists()


    @mock.patch(f"{snake_wrapper.__name__}.os.fork")
    def test_run_commandline_set_test_run(self, mock_fork):
        """Start a run in test mode, overriding config."""
        mock_fork.return_value = False
        expected_indir = (pathlib.Path(__file__).parent / "data" / "monitor_run" / "miniondir"
                          / "test1" / "rawdata" / "test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = str(pathlib.Path(__file__).parent / "data"
                         / "utilities_test" / "test_16s_config_no_debug.yaml")
        expected_outdir = (pathlib.Path(self.active_config["paths"]["output_base_path"]).relative_to(pathlib.Path(__file__).parent.parent)
                          / "NANO_Amplicon_Y20990101_RUN0001_XYZ-16S")
        assert not expected_outdir.exists()
        nomad_command = ["nomad", "job", "dispatch",
                         "-meta", f"indir={expected_indir}",
                         "-meta", f"outdir={expected_outdir}",
                         "-meta", f"runsheet={runsheet}",
                         "16s-snake-emu-staging", str(configfile)]
        expected_command = ["echo", f'"{" ".join(nomad_command)}"']
        test_args = ["--rundir", str(expected_indir), "--runsheet", str(runsheet),
                     "--workflow_config_file", str(configfile), "--test_run", "--dry_run"]
        test_command = snake_wrapper.run_pipeline(test_args)
        assert expected_command == test_command.args
        assert expected_outdir.exists()

    @mock.patch.dict(f"{snake_wrapper.__name__}.WORKFLOW_CONFIG", active_config,
                     clear = True)
    @mock.patch(f"{snake_wrapper.__name__}.DEFAULT_CONFIG_FILE", active_config_file)
    @mock.patch(f"{snake_wrapper.__name__}.os.fork")
    def test_fail_existing_dir(self, mock_fork):
        """Fail if the specified output directory already exists."""
        mock_fork.return_value = False
        expected_indir = (pathlib.Path(__file__).parent / "data" / "monitor_run" / "miniondir"
                          / "test1" / "rawdata" / "test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        expected_outdir = (pathlib.Path(__file__).parent / "data"
                           / "utilities_test"
                           / "test_existing_output")
        assert expected_outdir.exists()
        test_args = ["--rundir", str(expected_indir), "--runsheet", str(runsheet),
                     "--outdir", str(expected_outdir), "--dry_run"]
        error_msg = "The desired output directory already exists."
        with pytest.raises(FileExistsError, match=re.escape(error_msg)):
            snake_wrapper.run_pipeline(test_args)

    @mock.patch.dict(f"{snake_wrapper.__name__}.WORKFLOW_CONFIG", active_config,
                     clear = True)
    @mock.patch(f"{snake_wrapper.__name__}.DEFAULT_CONFIG_FILE", active_config_file)
    @mock.patch(f"{snake_wrapper.__name__}.os.fork")
    def test_continue_pipeline(self, mock_fork):
        """Continue the pipeline if specified in commandline mode."""
        mock_fork.return_value = False
        expected_indir = (pathlib.Path(__file__).parent / "data" / "monitor_run" / "miniondir"
                          / "test1" / "rawdata" / "test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        expected_outdir = (pathlib.Path(__file__).parent / "data"
                           / "utilities_test"
                           / "test_existing_output")
        configfile = self.active_config_file
        assert expected_outdir.exists()
        nomad_command = ["nomad", "job", "dispatch",
                         "-meta", f"indir={expected_indir}",
                         "-meta", f"outdir={expected_outdir}",
                         "-meta", f"runsheet={runsheet}",
                         "16s-snake-emu-staging", str(configfile)]
        expected_command = ["echo", f'"{" ".join(nomad_command)}"']
        test_args = ["--rundir", str(expected_indir), "--runsheet", str(runsheet),
                     "--outdir", str(expected_outdir),
                    "--continue_pipeline", "--dry_run"]
        test_command = snake_wrapper.run_pipeline(test_args)
        assert expected_command == test_command.args

    @mock.patch.dict(f"{snake_wrapper.__name__}.WORKFLOW_CONFIG", active_config,
                     clear = True)
    @mock.patch(f"{snake_wrapper.__name__}.DEFAULT_CONFIG_FILE", active_config_file)
    @mock.patch(f"{snake_wrapper.__name__}.monitor_run.start_on_file_found")
    @mock.patch(f"{snake_wrapper.__name__}.os.fork")
    def test_check_wait_time(self, mock_fork, mock_start):
        """Ensure the pipeline handles wait times over 24 hours correctly."""
        mock_fork.return_value = False
        # check we're waiting correctly
        watch_hours = 24
        fudge_hours = 1
        watch_seconds = (watch_hours + fudge_hours) * 3600

        expected_indir = (pathlib.Path(__file__).parent / "data" / "monitor_run" / "miniondir"
                          / "test1" / "rawdata" / "test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = self.active_config_file
        expected_outdir = pathlib.Path(self.active_config["paths"]["output_base_path"]) \
                          / "NANO_Amplicon_Y20990101_RUN0001_XYZ-16S"
        assert not expected_outdir.exists()
        nomad_command = ["nomad", "job", "dispatch",
                         "-meta", f"indir={expected_indir}",
                         "-meta", f"outdir={expected_outdir}",
                         "-meta", f"runsheet={runsheet}",
                         "16s-snake-emu-staging", str(configfile)]
        expected_command = ["echo", f'"{" ".join(nomad_command)}"']
        test_args = ["--rundir", str(expected_indir), "--runsheet", str(runsheet),
                     "--run_time", str(watch_hours), "--dry_run"]
        test_command = snake_wrapper.run_pipeline(test_args)
        assert expected_outdir.exists()
        mock_start.assert_called_with(mock.ANY, # run params - TODO: should we double-check we're calling with the right ones?
                                      "final_summary*.txt",
                                      dry_run=True,
                                      watch_timeout = watch_seconds,
                                      watch_interval = mock.ANY,
                                      log_interval=3600)

    @mock.patch(f"{snake_wrapper.__name__}.os.fork")
    def test_use_lis_from_file(self, mock_fork):
        """Use a LIS report from a file."""
        mock_fork.return_value = False
        expected_indir = (pathlib.Path(__file__).parent / "data" / "monitor_run" / "miniondir"
                          / "test1" / "rawdata" / "test_subdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = str(pathlib.Path(__file__).parent / "data"
                         / "utilities_test" / "test_16s_config_use_lis.yaml")
        expected_outdir = (pathlib.Path(self.active_config["paths"]["output_base_path"]).relative_to(pathlib.Path(__file__).parent.parent)
                          / "NANO_Amplicon_Y20990101_RUN0001_XYZ-16S")
        assert not expected_outdir.exists()
        nomad_command = ["nomad", "job", "dispatch",
                         "-meta", f"indir={expected_indir}",
                         "-meta", f"outdir={expected_outdir}",
                         "-meta", f"runsheet={runsheet}",
                         "16s-snake-emu-staging", str(configfile)]
        expected_command = ["echo", f'"{" ".join(nomad_command)}"']
        test_args = ["--rundir", str(expected_indir), "--runsheet", str(runsheet),
                     "--workflow_config_file", str(configfile), "--dry_run"]
        log_msg = "The runsheet is correct."
        with self._caplog.at_level(logging.INFO, logger="check_runsheet"):
            test_command = snake_wrapper.run_pipeline(test_args)
            assert ("check_runsheet", logging.INFO, log_msg) in self._caplog.record_tuples
        assert expected_command == test_command.args
        assert expected_outdir.exists()
