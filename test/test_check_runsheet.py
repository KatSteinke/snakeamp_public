import logging
import pathlib
import re
import unittest
from unittest import mock

import pandas as pd
import pytest

import check_runsheet


class TestCheckSinglePrefix(unittest.TestCase):
    test_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_translate_runsheet.xlsx"
    fake_mads = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "fake_mads_data.csv"
    sheet_data = pd.read_excel(test_runsheet, usecols = "A:B", skiprows = 3,
                               dtype = {"Prøvenummer": str})
    sheet_data = sheet_data.dropna()

    lab_info_data = pd.read_csv(fake_mads, encoding="latin1", dtype={"afsendt": str, "cprnr.": str,
                                                                     "modtaget": str})
    test_config = {"sample_number_settings": {"sample_number_format":
                                                  '([BDFPT]|[1357]0)([0-9]{8}|[0-9]{6})',
                                              "sample_numbers_in": "number",
                                              "sample_numbers_out": "letter",
                                              "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                              "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                              "number_to_letter": {"70": "P",
                                                                   "30": "B",
                                                                   "10": "D",
                                                                   "11": "F",
                                                                   "50": "T"},
                                              "negative_control": '',
                                              "positive_control": {}},
                   "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                   "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                   }

    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    def test_prefix_not_in_sheet(self):
        with pytest.raises(ValueError, match="No samples with prefix X found in runsheet."):
            check_runsheet.check_by_prefix(self.sheet_data, self.lab_info_data, "X", "P",
                                           active_config = self.test_config)

    def test_sample_not_in_report(self):
        fail_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test"\
                        / "test_notinmads_runsheet.xlsx"
        fail_data = pd.read_excel(fail_runsheet, usecols="A:C", skiprows=3,
                                  dtype={"Prøvenummer": str, "Barkode": str})
        fail_data = fail_data.dropna()
        error_msg = "Samples ['1121400000', '1121410000'] were not found in MADS report. " \
                    "Please check that sample numbers are correct."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_by_prefix(fail_data, self.lab_info_data, "11", "F",
                                           active_config = self.test_config)

    def test_catch_lis_duplicates(self):
        lis_with_duplicates = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" /\
                              "fake_mads_duplicated.csv"
        lab_info_data = pd.read_csv(lis_with_duplicates, encoding = "latin1",
                                    dtype = {"afsendt": str, "cprnr.": str,
                                             "modtaget": str})
        error_msg = "MADS report contains duplicated sample numbers. " \
                    "This likely means the report covers multiple years. " \
                    "Get a new MADS report with the correct start date."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_by_prefix(self.sheet_data, lab_info_data, "11", "F",
                                           active_config = self.test_config)

    def test_success(self):
        test_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_translate_runsheet.xlsx"
        sheet_data = pd.read_excel(test_runsheet, usecols = "A:B", skiprows = 3,
                                   dtype = {"Prøvenummer": str})
        sheet_data = sheet_data.dropna()
        success_msg = "All samples with prefix 11 found in LIS."
        with self._caplog.at_level(logging.DEBUG, logger = "check_runsheet"):
            check_runsheet.check_by_prefix(sheet_data, self.lab_info_data, "11", "F",
                                           active_config = self.test_config)

        assert ("check_runsheet", logging.DEBUG, success_msg) in self._caplog.record_tuples

    def test_success_new_format(self):
        fake_mads = pathlib.Path(
            __file__).parent / "data" / "sample_sheet_test" / "fake_mads_new_format.csv"

        lab_info_data = pd.read_csv(fake_mads, encoding = "latin1",
                                    dtype = {"afsendt": str, "cprnr.": str,
                                             "modtaget": str})
        test_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_translate_runsheet.xlsx"
        sheet_data = pd.read_excel(test_runsheet, usecols = "A:B", skiprows = 3,
                                   dtype = {"Prøvenummer": str})
        sheet_data = sheet_data.dropna()
        success_msg = "All samples with prefix 11 found in LIS."
        with self._caplog.at_level(logging.DEBUG, logger = "check_runsheet"):
            check_runsheet.check_by_prefix(sheet_data, lab_info_data, "11", "F",
                                           active_config = self.test_config)

        assert ("check_runsheet", logging.DEBUG, success_msg) in self._caplog.record_tuples

    def test_success_controls(self):
        """Ensure comparison against controls is performed"""
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      r'([BDFPT]|[1357]0)([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "11": "F",
                                                                       "50": "T"},
                                                  
                                                  "negative_control": 'NegK',
                                                  "positive_control": {"PosK": "Placeholderia"}},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        test_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_translate_runsheet.xlsx"
        sheet_data = pd.read_excel(test_runsheet, usecols = "A:B", skiprows = 3,
                                   dtype = {"Prøvenummer": str})
        sheet_data = sheet_data.dropna()
        success_msg = "All samples with prefix 11 found in LIS."
        with self._caplog.at_level(logging.DEBUG, logger = "check_runsheet"):
            check_runsheet.check_by_prefix(sheet_data, self.lab_info_data, "11", "F",
                                           active_config = test_config)

        assert ("check_runsheet", logging.DEBUG, success_msg) in self._caplog.record_tuples


class TestCheckRunsheetFormat(unittest.TestCase):
    test_config = {"sample_number_settings": {"sample_number_format":
                                                  '([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})',
                                              "sample_numbers_in": "number",
                                              "sample_numbers_out": "letter",
                                              # TODO: how to handle splicing in year?
                                              "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})?(?P<sample_number>\d{6})',
                                              "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                              "number_to_letter": {"70": "P",
                                                                   "30": "B",
                                                                   "11": "F",
                                                                   "10": "D",
                                                                   "50": "T"},
                                              
                                              "negative_control": '',
                                              "positive_control": {}},
                   "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                   "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                   }
    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    def test_multiple_fails(self):
        fail_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_notinmads_runsheet.xlsx"
        fake_mads = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "fake_mads_data.csv"
        sheet_data = pd.read_excel(fail_runsheet, usecols="A:C", skiprows=3,
                                   dtype = {"Prøvenummer": str})
        sheet_data = sheet_data.dropna()
        error_msg = "Samples ['1121400000', '1121410000'] were not found in MADS report. " \
                    "Please check that sample numbers are correct."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_against_lis(sheet_data, fake_mads,
                                             active_config = self.test_config)

    def test_success(self):
        test_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_translate_runsheet.xlsx"
        fake_mads = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "fake_mads_data.csv"
        sheet_data = pd.read_excel(test_runsheet, usecols = "A:C", skiprows = 3,
                                   dtype = {"Prøvenummer": str})
        sheet_data = sheet_data.dropna()
        success_msg = "The runsheet is correct."
        with self._caplog.at_level(logging.INFO, logger = "check_runsheet"):
            check_runsheet.check_against_lis(sheet_data, fake_mads, active_config = self.test_config)

        assert ("check_runsheet", logging.INFO, success_msg) in self._caplog.record_tuples

    def test_success_new_format(self):
        """Successfully handle a LIS report with new columns."""
        test_runsheet = (pathlib.Path(__file__).parent / "data" / "sample_sheet_test"
                         / "test_translate_runsheet.xlsx")
        fake_mads = (pathlib.Path(__file__).parent / "data" / "sample_sheet_test"
                     / "fake_mads_new_format.csv")
        sheet_data = pd.read_excel(test_runsheet, usecols = "A:C", skiprows = 3,
                                   dtype = {"Prøvenummer": str})
        sheet_data = sheet_data.dropna()
        success_msg = "The runsheet is correct."
        with self._caplog.at_level(logging.INFO, logger = "check_runsheet"):
            check_runsheet.check_against_lis(sheet_data, fake_mads,
                                             active_config = self.test_config)

        assert ("check_runsheet", logging.INFO, success_msg) in self._caplog.record_tuples

    def test_success_controls(self):
        test_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_translate_runsheet.xlsx"
        fake_mads = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "fake_mads_data.csv"
        sheet_data = pd.read_excel(test_runsheet, usecols = "A:C", skiprows = 3,
                                   dtype = {"Prøvenummer": str})
        sheet_data = sheet_data.dropna()
        success_msg = "The runsheet is correct."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})?(?P<sample_number>\d{6})',
                                                  "format_in_lis": r'(?P<sample_type>[BDPFT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "11": "F",
                                                                       "10": "D",
                                                                       "50": "T"},
                                                  
                                                  "negative_control": 'NegK',
                                                  "positive_control": {"PosK": "Placeholderia"}},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        with self._caplog.at_level(logging.INFO, logger = "check_runsheet"):
            check_runsheet.check_against_lis(sheet_data, fake_mads, active_config = test_config)
            assert ("check_runsheet", logging.INFO, success_msg) in self._caplog.record_tuples

    def test_success_rearrange(self):
        test_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_translate_runsheet_rearrange.xlsx"
        fake_mads = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "fake_mads_data_move_year.csv"
        sheet_data = pd.read_excel(test_runsheet, usecols = "A:C", skiprows = 3,
                                   dtype = {"Prøvenummer": str})
        sheet_data = sheet_data.dropna()
        success_msg = "The runsheet is correct."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_number>\d{6})(?P<sample_year>\d{2})',
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "11": "F",
                                                                       "10": "D",
                                                                       "50": "T"},
                                                  
                                                  "negative_control": '',
                                                  "positive_control": {}},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        with self._caplog.at_level(logging.INFO, logger = "check_runsheet"):
            check_runsheet.check_against_lis(sheet_data, fake_mads, active_config = test_config)
            assert ("check_runsheet", logging.INFO, success_msg) in self._caplog.record_tuples

    def test_success_drop_component(self):
        """Alert the user when sample number format in LIS report contains fewer components than
        sample number format in sheet."""
        test_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_translate_runsheet_isolate.xlsx"
        fake_mads = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "fake_mads_data.csv"
        sheet_data = pd.read_excel(test_runsheet, usecols = "A:C", skiprows = 3,
                                   dtype = {"Prøvenummer": str})
        sheet_data = sheet_data.dropna()
        dropped_component_msg = ("Comparing only"
                                 " ['sample_type', 'sample_year', 'sample_number'] to LIS report. "
                                 "Cannot check if ['bact_number'] component(s) are correct.")
        success_msg = "The runsheet is correct."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      r'([BFDPT]|[1357]0|11)([0-9]{8}|[0-9]{6})(-\d)?',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                  "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "11": "F",
                                                                       "10": "D",
                                                                       "50": "T"},
                                                  
                                                  "negative_control": '',
                                                  "positive_control": {}},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        with self._caplog.at_level(logging.INFO, logger = "check_runsheet"):
            check_runsheet.check_against_lis(sheet_data, fake_mads, active_config = test_config)
            assert ("check_runsheet", logging.INFO, success_msg) in self._caplog.record_tuples
            assert ("check_runsheet", logging.INFO, dropped_component_msg) in self._caplog.record_tuples


class TestCheckSampleNumbers(unittest.TestCase):
    test_config = {"sample_number_settings": {"sample_number_format":
                                                  '([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})',
                                              "sample_numbers_in": "number",
                                              "sample_numbers_out": "letter",
                                              "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_number>\d{6})',
                                              "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                              "number_to_letter": {"70": "P",
                                                                   "30": "B",
                                                                   "10": "D",
                                                                   "11": "F",
                                                                   "50": "T"},
                                              
                                              "negative_control": '',
                                              "positive_control": {}},
                   "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                   "barcode_prefix": "RB"
                   # barcode prefix as letter (for transferring original fastqs by barcode)
                   }
    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    def test_fail_ids(self):
        id_fail_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-id-fail.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:\n" \
                    "Sample IDs ['123'] are not valid. " \
                    "Sample IDs must start with 70 or 30 or 10 or 11 or 50 followed by eight numbers" \
                    " (six if leaving out year). " \
                    "Please correct sample IDs in runsheet."
        sheet_data = pd.read_excel(id_fail_sheet, usecols="A:B", skiprows=3,
                                   dtype={"Prøvenummer": str})
        sheet_data = sheet_data.dropna()
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data, active_config = self.test_config)

    def test_fail_ids_letters(self):
        id_fail_sheet = pathlib.Path(__file__).parent /"data"/ "utilities_test" \
                        / "runsheet-letters.xlsx"
        sheet_data = pd.read_excel(id_fail_sheet, usecols="A:B", skiprows=3,
                                   dtype={"Prøvenummer": str})
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nSample IDs ['1199123456'] are not valid." \
                    " Sample IDs must start with P or B or D or F or T followed by eight numbers" \
                    " (six if leaving out year). " \
                    "Please correct sample IDs in runsheet."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '[BDFPT]([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "letter",
                                                  "sample_numbers_out": "letter",
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "11": "F",
                                                                       "50": "T"},
                                                  
                                                  "negative_control": '',
                                                  "positive_control": {}},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data, active_config = test_config)

    def test_id_fail_negk(self):
        fail_id_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                        / "runsheet-id-fail-negk.xlsx"
        sheet_data = pd.read_excel(fail_id_sheet, usecols = "A:B", skiprows = 3,
                                   dtype = {"Prøvenummer": str})
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nSample IDs ['123'] are not valid." \
                    " Sample IDs must start with 70 or 30 or 10 or 11 or 50 followed by eight numbers" \
                    " (six if leaving out year). " \
                    "Negative controls must be given in the format NegK. " \
                    "Please correct sample IDs in runsheet."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "11": "F",
                                                                       "50": "T"},
                                                  
                                                  "negative_control": 'NegK',
                                                  "positive_control": {}},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data, active_config = test_config)

    def test_fail_no_ids(self):
        no_id_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-no-id.xlsx"
        sheet_data = pd.read_excel(no_id_sheet, usecols = "A:B", skiprows = 3,
                                   dtype = {"Prøvenummer": str})
        sheet_data = sheet_data.dropna()
        error_msg = "The following issue(s) were detected with the runsheet:\n" \
                    "No sample IDs found."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data, active_config = self.test_config)

    def test_fail_no_positive_control(self):
        no_positive_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-no-posk.xlsx"
        sheet_data = pd.read_excel(no_positive_sheet, usecols="A:B", skiprows=3,
                                   dtype={"Prøvenummer": str})
        sheet_data = sheet_data.dropna()
        error_msg = "The following issue(s) were detected with the runsheet:\n" \
                    "No positive controls given in runsheet."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "11": "F",
                                                                       "50": "T"},
                                                  
                                                  "negative_control": '',
                                                  "positive_control": {"PosK": "Placeholderia"}},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data, active_config = test_config)

    def test_fail_no_negative_control(self):
        no_negative_sheet = pathlib.Path(__file__).parent / "data" /"utilities_test" / "runsheet-no-negk.xlsx"
        sheet_data = pd.read_excel(no_negative_sheet, usecols="A:B", skiprows=3,
                                   dtype={"Prøvenummer": str})
        sheet_data = sheet_data.dropna()
        error_msg = "No negative controls given in runsheet."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})',
                                                  "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "11": "F",
                                                                       "50": "T"},
                                                  
                                                  "negative_control": 'NegK',
                                                  "positive_control": {}},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data, active_config = test_config)

    def test_warn_duplicated_ids(self):
        duplicated_id_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-id-duplication.xlsx"
        sheet_data = pd.read_excel(duplicated_id_sheet, usecols="A:B", skiprows=3,
                                   dtype={"Prøvenummer": str})
        duplicated_warning = "Sample number(s) ['1123456789'] are duplicated." \
                             " If you are sure you want to sequence the same sample twice, " \
                             "you can ignore this warning."
        with self._caplog.at_level(logging.WARNING, logger = "check_runsheet"):
            check_runsheet.check_sheet_format(sheet_data, active_config = self.test_config)
            assert ("check_runsheet", logging.WARNING, duplicated_warning) in self._caplog.record_tuples

    def test_fail_barcodes(self):
        barcode_fail_sheet = pathlib.Path(__file__).parent /"data" /"utilities_test" / "runsheet-barcode-fail.xlsx"
        sheet_data = pd.read_excel(barcode_fail_sheet, usecols = "A:C", skiprows = 3,
                                   dtype = {"Prøvenummer": str})
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nBarcodes ['RB3'] are not valid barcodes. " \
                    "Barcodes must consist of RB + a number between 01 and 96."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data, check_barcodes = True,
                                              active_config = self.test_config)

    def test_fail_no_barcodes(self):
        no_barcode_sheet = pathlib.Path(__file__).parent / "data" /"utilities_test" / "runsheet-no-barcode.xlsx"
        sheet_data = pd.read_excel(no_barcode_sheet, usecols = "A:C", skiprows = 3,
                                   dtype = {"Prøvenummer": str})
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nNo barcodes found."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data, check_barcodes = True,
                                              active_config = self.test_config)

    def test_fail_more_barcodes(self):
        more_barcodes_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-more-barcodes.xlsx"
        sheet_data = pd.read_excel(more_barcodes_sheet, usecols = "A:C", skiprows = 3,
                                   dtype = {"Prøvenummer": str})
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nAmount of sample IDs and barcodes don't match. " \
                    "There are 2 sample IDs but 3 barcodes."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data, check_barcodes = True,
                                              active_config = self.test_config)

    def test_fail_duplicated_barcodes(self):
        duplicated_barcodes_sheet = pathlib.Path(
            __file__).parent / "data" / "utilities_test" / "runsheet-barcode-duplication.xlsx"
        sheet_data = pd.read_excel(duplicated_barcodes_sheet, usecols = "A:C", skiprows = 3,
                                   dtype = {"Prøvenummer": str})
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nBarcode(s) ['RB02'] are duplicated."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data, check_barcodes = True,
                                              active_config = self.test_config)


class TestCheckRunsheet(unittest.TestCase):
    test_config = {"sample_number_settings": {"sample_number_format":
                                                  '([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})',
                                              "sample_numbers_in": "number",
                                              "sample_numbers_out": "letter",
                                              "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_number>\d{6})',
                                              "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                              "number_to_letter": {"70": "P",
                                                                   "30": "B",
                                                                   "10": "D",
                                                                   "11": "F",
                                                                   "50": "T"},
                                              
                                              "negative_control": '',
                                              "positive_control": {}},
                   "lab_info_system": {"use_lis_features": False},
                   "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                   "barcode_prefix": "RB"
                   # barcode prefix as letter (for transferring original fastqs by barcode)
                   }

    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    def test_fail_ids(self):
        id_fail_sheet = pathlib.Path(
            __file__).parent / "data" / "utilities_test" / "runsheet-id-fail.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:\n" \
                    "Sample IDs ['123'] are not valid. " \
                    "Sample IDs must start with 70 or 30 or 10 or 11 or 50 followed by eight numbers" \
                    " (six if leaving out year). " \
                    "Please correct sample IDs in runsheet."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            check_runsheet.check_runsheet(id_fail_sheet, active_config = self.test_config)

    def test_fail_ids_letters(self):
        id_fail_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                        / "runsheet-letters.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nSample IDs ['1199123456'] are not valid." \
                    " Sample IDs must start with P or B or D or F or T followed by eight numbers" \
                    " (six if leaving out year). " \
                    "Please correct sample IDs in runsheet."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '[BDFPT]([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "letter",
                                                  "sample_numbers_out": "letter",
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "11": "F",
                                                                       "50": "T"},
                                                  
                                                  "negative_control": '',
                                                  "positive_control": {}},
                       "lab_info_system": {"use_lis_features": False},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"
                       # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            check_runsheet.check_runsheet(id_fail_sheet, active_config = test_config)

    def test_id_fail_negk(self):
        fail_id_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                        / "runsheet-id-fail-negk.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nSample IDs ['123'] are not valid." \
                    " Sample IDs must start with 70 or 30 or 10 or 11 or 50 followed by eight numbers" \
                    " (six if leaving out year). " \
                    "Negative controls must be given in the format NegK. " \
                    "Please correct sample IDs in runsheet."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "11": "F",
                                                                       "50": "T"},
                                                  
                                                  "negative_control": 'NegK',
                                                  "positive_control": {}},
                       "lab_info_system": {"use_lis_features": False},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"
                       # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            check_runsheet.check_runsheet(fail_id_sheet, active_config = test_config)

    def test_fail_no_ids(self):
        no_id_sheet = pathlib.Path(
            __file__).parent / "data" / "utilities_test" / "runsheet-no-id.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:\n" \
                    "No sample IDs found."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            check_runsheet.check_runsheet(no_id_sheet, active_config = self.test_config)

    def test_fail_no_positive_control(self):
        no_positive_sheet = pathlib.Path(
            __file__).parent / "data" / "utilities_test" / "runsheet-no-posk.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:\n" \
                    "No positive controls given in runsheet."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "11": "F",
                                                                       "50": "T"},
                                                  
                                                  "negative_control": '',
                                                  "positive_control": {"PosK": "Placeholderia"}},
                       "lab_info_system": {"use_lis_features": False},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"
                       # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            check_runsheet.check_runsheet(no_positive_sheet, active_config = test_config)

    def test_fail_no_negative_control(self):
        no_negative_sheet = pathlib.Path(
            __file__).parent / "data" / "utilities_test" / "runsheet-no-negk.xlsx"
        error_msg = "No negative controls given in runsheet."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})',
                                                  "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "11": "F",
                                                                       "50": "T"},
                                                  
                                                  "negative_control": 'NegK',
                                                  "positive_control": {}},
                       "lab_info_system": {"use_lis_features": False},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"
                       # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            check_runsheet.check_runsheet(no_negative_sheet, active_config = test_config)

    def test_fail_barcodes(self):
        barcode_fail_sheet = pathlib.Path(__file__).parent /"data" /"utilities_test" / "runsheet-barcode-fail.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nBarcodes ['RB3'] are not valid barcodes. " \
                    "Barcodes must consist of RB + a number between 01 and 96."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_runsheet(barcode_fail_sheet, check_barcodes = True,
                                          active_config = self.test_config)

    def test_fail_no_barcodes(self):
        no_barcode_sheet = pathlib.Path(__file__).parent / "data" /"utilities_test" / "runsheet-no-barcode.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nNo barcodes found."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_runsheet(no_barcode_sheet, check_barcodes = True,
                                          active_config = self.test_config)

    def test_fail_more_barcodes(self):
        more_barcodes_sheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                               / "runsheet-more-barcodes.xlsx")

        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nAmount of sample IDs and barcodes don't match. " \
                    "There are 2 sample IDs but 3 barcodes."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_runsheet(more_barcodes_sheet, check_barcodes = True,
                                          active_config = self.test_config)

    def test_fail_duplicated_barcodes(self):
        duplicated_barcodes_sheet = pathlib.Path(
            __file__).parent / "data" / "utilities_test" / "runsheet-barcode-duplication.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nBarcode(s) ['RB02'] are duplicated."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            check_runsheet.check_runsheet(duplicated_barcodes_sheet, check_barcodes = True,
                                          active_config = self.test_config)

    def test_missing_from_lis(self):
        """Fail if the LIS report should be checked and the sample is missing from it."""
        runsheet = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                   / "test_nanopore_runsheet.xlsx"
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      r'([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})(-\d)?',
                                                  "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                  "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "sample_numbers_in": "letter",
                                                  "sample_numbers_out": "letter",
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "11": "F",
                                                                       "50": "T"},
                                                  
                                                  "negative_control": 'NegK[a-zA-Z0-9_-]*',
                                                  "positive_control": {}},
                       "lab_info_system": {"use_lis_features": True,
                                           "lis_report": str(pathlib.Path(__file__).parent
                                                             / "data" / "sample_sheet_test"
                                                             / "fake_mads_data.csv")},
                       "barcode_format": "NB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "NB"
                       # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        error_msg = "Samples ['F99123456-1', 'F99123456-2'] were not found in MADS report. " \
                    "Please check that sample numbers are correct."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            check_runsheet.check_runsheet(runsheet, active_config = test_config)

    def test_success_no_lis(self):
        """Successfully check the runsheet without using a LIS report."""
        runsheet = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                   / "test_nanopore_runsheet.xlsx"
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      r'([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})(-\d)?',
                                                  "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                  "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "sample_numbers_in": "letter",
                                                  "sample_numbers_out": "letter",
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "11": "F",
                                                                       "50": "T"},
                                                  
                                                  "negative_control": 'NegK[a-zA-Z0-9_-]*',
                                                  "positive_control": {}},
                       "lab_info_system": {"use_lis_features": False,
                                           "lis_report": str(pathlib.Path(__file__).parent
                                                             / "data" / "sample_sheet_test"
                                                             / "fake_mads_data.csv")},
                       "barcode_format": "NB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "NB"
                       # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        loading_sheet_msg = "Loading runsheet..."
        check_sheet_msg = "Checking runsheet format...."
        check_lis_msg = "Comparing to samples in MADS......"
        with self._caplog.at_level(logging.INFO, logger = "check_runsheet"):
            runsheet_pass = check_runsheet.check_runsheet(runsheet, active_config = test_config)

        assert ("check_runsheet", logging.INFO, loading_sheet_msg) in self._caplog.record_tuples
        assert ("check_runsheet", logging.INFO, check_sheet_msg) in self._caplog.record_tuples
        assert ("check_runsheet", logging.INFO, check_lis_msg) not in self._caplog.record_tuples
        assert runsheet_pass

    def test_success_use_lis(self):
        """Successfully check the runsheet while using a LIS report."""
        runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" \
                   / "test_translate_runsheet.xlsx"
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      r'([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})(-\d)?',
                                                  "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                  "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "11": "F",
                                                                       "50": "T"},
                                                  "date_settings":
                                                      {"splice_in_date": True,
                                                       "length_without_date": 8,
                                                       "splice_after": 2},
                                                  "negative_control": '',
                                                  "positive_control": {}},
                       "lab_info_system": {"use_lis_features": True,
                                           "lis_report": str(pathlib.Path(__file__).parent
                                                             / "data" / "sample_sheet_test"
                                                             / "fake_mads_data.csv")},
                       "barcode_format": "NB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "NB"
                       # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        loading_sheet_msg = "Loading runsheet..."
        check_sheet_msg = "Checking runsheet format...."
        check_lis_msg = "Comparing to samples in MADS......"
        with self._caplog.at_level(logging.INFO, logger = "check_runsheet"):
            runsheet_pass = check_runsheet.check_runsheet(runsheet, active_config = test_config)

        assert ("check_runsheet", logging.INFO, loading_sheet_msg) in self._caplog.record_tuples
        assert ("check_runsheet", logging.INFO, check_sheet_msg) in self._caplog.record_tuples
        assert ("check_runsheet", logging.INFO, check_lis_msg) in self._caplog.record_tuples
        assert runsheet_pass

    def test_success_use_new_lis(self):
        """Successfully check the runsheet while using a LIS report with additional columns."""
        runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" \
                   / "test_translate_runsheet.xlsx"
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      r'([BDFPT]|[1357]0|11)([0-9]{8}|[0-9]{6})(-\d)?',
                                                  "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                  "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "11": "F",
                                                                       "50": "T"},
                                                  "negative_control": '',
                                                  "positive_control": {}},
                       "lab_info_system": {"use_lis_features": True,
                                           "lis_report": str(pathlib.Path(__file__).parent
                                                             / "data" / "sample_sheet_test"
                                                             / "fake_mads_new_format.csv")},
                       "barcode_format": "NB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "NB"
                       # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        loading_sheet_msg = "Loading runsheet..."
        check_sheet_msg = "Checking runsheet format...."
        check_lis_msg = "Comparing to samples in MADS......"
        with self._caplog.at_level(logging.INFO, logger = "check_runsheet"):
            runsheet_pass = check_runsheet.check_runsheet(runsheet, active_config = test_config)

        assert ("check_runsheet", logging.INFO, loading_sheet_msg) in self._caplog.record_tuples
        assert ("check_runsheet", logging.INFO, check_sheet_msg) in self._caplog.record_tuples
        assert ("check_runsheet", logging.INFO, check_lis_msg) in self._caplog.record_tuples
        assert runsheet_pass


class TestRunCheck(unittest.TestCase):
    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    @mock.patch("builtins.input", side_effect = [str(pathlib.Path(__file__).parent / "data" / "utilities_test" \
                   / "test_nanopore_runsheet.xlsx")])
    def test_success_classic(self, mock_input):
        """Successfully check a good runsheet that has been entered by the user."""
        check_msg = "###Runsheet check"
        with (pytest.raises(SystemExit, match="0"),
            self._caplog.at_level(level="INFO", logger="check_runsheet")):
            check_runsheet.run_check(["--workflow_config_file", str(pathlib.Path(__file__).parent / "data"/"utilities_test"/"test_16s_config_no_debug.yaml")])
        assert ("check_runsheet", logging.INFO, check_msg) in self._caplog.record_tuples


    def test_success_commandline(self):
        """Successfully check a good runsheet that has been supplied as a commandline argument."""
        with (pytest.raises(SystemExit, match = "0")):
            check_runsheet.run_check(["--runsheet", str(pathlib.Path(__file__).parent / "data" / "utilities_test" \
                   / "test_nanopore_runsheet.xlsx"),
                                      "--workflow_config_file",
                                      str(pathlib.Path(__file__).parent / "data"/"utilities_test"/"test_16s_config_no_debug.yaml")])

    def test_catch_bad_sheet(self):
        """Successfully catch a broken runsheet."""
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nSample IDs ['1112345678', '1123456789', '123'] are not valid." \
                    " Sample IDs must start with P or B or D or F or T followed by eight numbers" \
                    " (six if leaving out year). " \
                    "Negative controls must be given in the format NegK[a-zA-Z0-9]*. " \
                    "Please correct sample IDs in runsheet."
        fail_runsheet = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                        / "runsheet-id-fail-negk.xlsx"
        with ((pytest.raises(SystemExit, match = "1")),
              self._caplog.at_level(level="INFO",logger="check_runsheet")):
            check_runsheet.run_check(["--runsheet",str(fail_runsheet),
                                      "--workflow_config_file",
                                      str(pathlib.Path(__file__).parent / "data"/"utilities_test"/"test_16s_config_no_debug.yaml")])
        assert ("check_runsheet", logging.ERROR, error_msg) in self._caplog.record_tuples

