import pathlib
import re
import unittest

from unittest import mock

import pandas as pd
import pytest

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


