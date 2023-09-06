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


class TestValidateRunsheet(unittest.TestCase):
    def test_fail_ids(self):
        id_fail_sheet = pathlib.Path(__file__).parent /"data"/ "utilities_test" / "runsheet-id-fail.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nSample IDs ['123'] are not valid." \
                    " Sample IDs must start with 70 or 30 or 10 or 50 followed by eight numbers" \
                    " (six if leaving out year). " \
                    "Please correct sample IDs in runsheet."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            snake_wrapper.validate_runsheet_format(id_fail_sheet)

    def test_fail_ids_letters(self):
        id_fail_sheet = pathlib.Path(__file__).parent /"data"/ "utilities_test" / "runsheet-letters.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nSample IDs ['1199123456'] are not valid." \
                    " Sample IDs must start with P or B or D or T followed by eight numbers" \
                    " (six if leaving out year). " \
                    "Please correct sample IDs in runsheet."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "letter",
                                                  "sample_numbers_out": "letter",
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "50": "T"},
                                                  "date_settings":
                                                      {"splice_in_date": False,
                                                       "length_without_date": 8,
                                                       "splice_after": 2},
                                                  "negative_control": '',
                                                  "positive_control": {}},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            snake_wrapper.validate_runsheet_format(id_fail_sheet, active_config = test_config)

    def test_fail_barcodes(self):
        barcode_fail_sheet = pathlib.Path(__file__).parent /"data" /"utilities_test" / "runsheet-barcode-fail.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nBarcodes ['RB3'] are not valid barcodes. " \
                    "Barcodes must consist of RB + a number between 01 and 96."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            snake_wrapper.validate_runsheet_format(barcode_fail_sheet)

    def test_fail_no_ids(self):
        no_id_sheet = pathlib.Path(__file__).parent / "data" /"utilities_test" / "runsheet-no-id.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nNo sample IDs found."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            snake_wrapper.validate_runsheet_format(no_id_sheet)

    def test_fail_no_barcodes(self):
        no_barcode_sheet = pathlib.Path(__file__).parent / "data" /"utilities_test" / "runsheet-no-barcode.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nNo barcodes found."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            snake_wrapper.validate_runsheet_format(no_barcode_sheet)

    def test_fail_more_barcodes(self):
        more_barcodes_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-more-barcodes.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nAmount of sample IDs and barcodes don't match. " \
                    "There are 2 sample IDs but 3 barcodes."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            snake_wrapper.validate_runsheet_format(more_barcodes_sheet)

    def test_fail_duplicated_barcodes(self):
        duplicated_barcodes_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-barcode-duplication.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nBarcode(s) ['RB02'] are duplicated."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            snake_wrapper.validate_runsheet_format(duplicated_barcodes_sheet)

    def test_fail_no_positive_control(self):
        no_positive_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-no-posk.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "No positive controls given in runsheet."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "50": "T"},
                                                  "date_settings":
                                                      {"splice_in_date": False,
                                                       "length_without_date": 8,
                                                       "splice_after": 2},
                                                  "negative_control": '',
                                                  "positive_control": {"PosK": "Placeholderia"}},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            snake_wrapper.validate_runsheet_format(no_positive_sheet, active_config = test_config)

    def test_fail_no_negative_control(self):
        no_negative_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-no-negk.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "No negative controls given in runsheet."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "50": "T"},
                                                  "date_settings":
                                                      {"splice_in_date": False,
                                                       "length_without_date": 8,
                                                       "splice_after": 2},
                                                  "negative_control": 'NegK',
                                                  "positive_control": {}},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            snake_wrapper.validate_runsheet_format(no_negative_sheet, active_config = test_config)
