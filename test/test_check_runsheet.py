import pathlib
import re
import unittest

import pandas as pd
import pytest

import check_runsheet


class TestCheckSinglePrefix(unittest.TestCase):
    test_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_translate_runsheet.xlsx"
    fake_mads = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "fake_mads_data.csv"
    sheet_data = pd.read_excel(test_runsheet, usecols = "A:B", skiprows = 3,
                               dtype = {"KMA nr": str})
    sheet_data = sheet_data.dropna()

    lab_info_data = pd.read_csv(fake_mads, encoding="latin1", dtype={"afsendt": str, "cprnr.": str,
                                                                     "modtaget": str})
    test_config = {"sample_number_settings": {"sample_number_format":
                                                  '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})',
                                              "sample_numbers_in": "number",
                                              "sample_numbers_out": "letter",
                                              "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_number>\d{6})',
                                              "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
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

    def test_prefix_not_in_sheet(self):
        with pytest.raises(ValueError, match="No samples with prefix X found in runsheet."):
            check_runsheet.check_by_prefix(self.sheet_data, self.lab_info_data, "X", "P",
                                           active_config = self.test_config)

    def test_sample_not_in_report(self):
        fail_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test"\
                        / "test_notinmads_runsheet.xlsx"
        fail_data = pd.read_excel(fail_runsheet, usecols="A:C", skiprows=3,
                                  dtype={"KMA nr": str, "Barkode NB": str})
        fail_data = fail_data.dropna()
        error_msg = "Samples ['11410000'] were not found in MADS report. " \
                    "Please check that sample numbers are correct."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_by_prefix(fail_data, self.lab_info_data, "30", "B",
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
            check_runsheet.check_by_prefix(self.sheet_data, lab_info_data, "30", "B",
                                           active_config = self.test_config)

    def test_success(self):
        test_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_translate_runsheet.xlsx"
        sheet_data = pd.read_excel(test_runsheet, usecols = "A:B", skiprows = 3,
                                   dtype = {"KMA nr": str})
        sheet_data = sheet_data.dropna()
        success_msg = "DEBUG:check_runsheet:All samples with prefix 30 found in LIS."
        with self.assertLogs("check_runsheet", level="DEBUG") as logged:
            check_runsheet.check_by_prefix(sheet_data, self.lab_info_data, "30", "B",
                                           active_config = self.test_config)
            assert success_msg in logged.output

    def test_success_controls(self):
        """Ensure comparison against controls is performed"""
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_number>\d{6})',
                                                  "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "50": "T"},
                                                  "date_settings":
                                                      {"splice_in_date": False,
                                                       "length_without_date": 8,
                                                       "splice_after": 2},
                                                  "negative_control": 'NegK',
                                                  "positive_control": {"PosK": "Placeholderia"}},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        test_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_translate_runsheet.xlsx"
        sheet_data = pd.read_excel(test_runsheet, usecols = "A:B", skiprows = 3,
                                   dtype = {"KMA nr": str})
        sheet_data = sheet_data.dropna()
        success_msg = "DEBUG:check_runsheet:All samples with prefix 30 found in LIS."
        with self.assertLogs("check_runsheet", level = "DEBUG") as logged:
            check_runsheet.check_by_prefix(sheet_data, self.lab_info_data, "30", "B",
                                           active_config = test_config)
            assert success_msg in logged.output


class TestCheckRunsheet(unittest.TestCase):
    test_config = {"sample_number_settings": {"sample_number_format":
                                                  '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})',
                                              "sample_numbers_in": "number",
                                              "sample_numbers_out": "letter",
                                              "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_number>\d{6})',
                                              "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
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

    def test_multiple_fails(self):
        fail_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_notinmads_runsheet.xlsx"
        fake_mads = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "fake_mads_data.csv"
        sheet_data = pd.read_excel(fail_runsheet, usecols="A:B", skiprows=3,
                                   dtype={"KMA nr": str})
        sheet_data = sheet_data.dropna()
        error_msg = "The following issues were encountered:\n" \
                    "Samples ['11400000'] were not found in MADS report. " \
                    "Please check that sample numbers are correct.\n" \
                    "Samples ['11410000'] were not found in MADS report. " \
                    "Please check that sample numbers are correct."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_against_lis(sheet_data, fake_mads,
                                             active_config = self.test_config)

    def test_catch_prefix_fail(self):
        fail_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_notinmads_runsheet.xlsx"
        fake_mads = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "fake_mads_data.csv"
        sheet_data = pd.read_excel(fail_runsheet, usecols = "A:B", skiprows = 3,
                                   dtype = {"KMA nr": str})
        sheet_data = sheet_data.dropna()
        error_msg = "The following issues were found with the runsheet:\n" \
                    "No valid prefixes found in runsheet."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            check_runsheet.check_against_lis(sheet_data, fake_mads)

    def test_success(self):
        test_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_translate_runsheet.xlsx"
        fake_mads = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "fake_mads_data.csv"
        sheet_data = pd.read_excel(test_runsheet, usecols = "A:B", skiprows = 3,
                                   dtype = {"KMA nr": str})
        sheet_data = sheet_data.dropna()
        success_msg = "INFO:check_runsheet:The runsheet is correct."
        with self.assertLogs("check_runsheet") as logged:
            check_runsheet.check_against_lis(sheet_data, fake_mads, active_config = self.test_config)
            assert success_msg in logged.output

    def test_success_controls(self):
        test_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_translate_runsheet.xlsx"
        fake_mads = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "fake_mads_data.csv"
        sheet_data = pd.read_excel(test_runsheet, usecols = "A:B", skiprows = 3,
                                   dtype = {"KMA nr": str})
        sheet_data = sheet_data.dropna()
        success_msg = "INFO:check_runsheet:The runsheet is correct."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_number>\d{6})',
                                                  "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "number_to_letter": {"70": "P",
                                                                       "30": "B",
                                                                       "10": "D",
                                                                       "50": "T"},
                                                  "date_settings":
                                                      {"splice_in_date": False,
                                                       "length_without_date": 8,
                                                       "splice_after": 2},
                                                  "negative_control": 'NegK',
                                                  "positive_control": {"PosK": "Placeholderia"}},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        with self.assertLogs("check_runsheet") as logged:
            check_runsheet.check_against_lis(sheet_data, fake_mads, active_config = test_config)
            assert success_msg in logged.output


class TestCheckSampleNumbers(unittest.TestCase):
    def test_fail_ids(self):
        id_fail_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-id-fail.xlsx"
        error_msg = "The following issue(s) were detected with the runsheet:\n" \
                    "Sample IDs ['123'] are not valid. " \
                    "Sample IDs must start with 70 or 30 or 10 or 50 followed by eight numbers" \
                    " (six if leaving out year). " \
                    "Please correct sample IDs in runsheet."
        sheet_data = pd.read_excel(id_fail_sheet, usecols="A:B", skiprows=3,
                                   dtype={"KMA nr": str})
        sheet_data = sheet_data.dropna()
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data)

    def test_fail_ids_letters(self):
        id_fail_sheet = pathlib.Path(__file__).parent /"data"/ "utilities_test" \
                        / "runsheet-letters.xlsx"
        sheet_data = pd.read_excel(id_fail_sheet, usecols="A:B", skiprows=3,
                                   dtype={"KMA nr": str})
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nSample IDs ['1199123456'] are not valid." \
                    " Sample IDs must start with P or B or D or T followed by eight numbers" \
                    " (six if leaving out year). " \
                    "Please correct sample IDs in runsheet."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '[BDPT]([0-9]{8}|[0-9]{6})',
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
            check_runsheet.check_sheet_format(sheet_data, active_config = test_config)

    def test_id_fail_negk(self):
        fail_id_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                        / "runsheet-id-fail-negk.xlsx"
        sheet_data = pd.read_excel(fail_id_sheet, usecols = "A:B", skiprows = 3,
                                   dtype = {"KMA nr": str})
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nSample IDs ['123'] are not valid." \
                    " Sample IDs must start with 70 or 30 or 10 or 50 followed by eight numbers" \
                    " (six if leaving out year). " \
                    "Negative controls must be given in the format NegK. " \
                    "Please correct sample IDs in runsheet."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "format_in_sheet": '(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "format_in_lis": '(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
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
            check_runsheet.check_sheet_format(sheet_data, active_config = test_config)

    def test_fail_no_ids(self):
        no_id_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-no-id.xlsx"
        sheet_data = pd.read_excel(no_id_sheet, usecols = "A:B", skiprows = 3,
                                   dtype = {"KMA nr": str})
        sheet_data = sheet_data.dropna()
        error_msg = "The following issue(s) were detected with the runsheet:\n" \
                    "No sample IDs found."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data)

    def test_fail_no_positive_control(self):
        no_positive_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-no-posk.xlsx"
        sheet_data = pd.read_excel(no_positive_sheet, usecols="A:B", skiprows=3,
                                   dtype={"KMA nr": str})
        sheet_data = sheet_data.dropna()
        error_msg = "The following issue(s) were detected with the runsheet:\n" \
                    "No positive controls given in runsheet."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})',
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
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
            check_runsheet.check_sheet_format(sheet_data, active_config = test_config)

    def test_fail_no_negative_control(self):
        no_negative_sheet = pathlib.Path(__file__).parent / "data" /"utilities_test" / "runsheet-no-negk.xlsx"
        sheet_data = pd.read_excel(no_negative_sheet, usecols="A:B", skiprows=3,
                                   dtype={"KMA nr": str})
        sheet_data = sheet_data.dropna()
        error_msg = "No negative controls given in runsheet."
        test_config = {"sample_number_settings": {"sample_number_format":
                                                      '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})',
                                                  "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
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
            check_runsheet.check_sheet_format(sheet_data, active_config = test_config)

    def test_warn_duplicated_ids(self):
        duplicated_id_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-id-duplication.xlsx"
        sheet_data = pd.read_excel(duplicated_id_sheet, usecols="A:B", skiprows=3,
                                   dtype={"KMA nr": str})
        with self.assertLogs("check_runsheet") as logged:
            check_runsheet.check_sheet_format(sheet_data)
            duplicated_warning = "WARNING:check_runsheet:Sample number(s) ['1123456789'] are duplicated." \
                                 " If you are sure you want to sequence the same sample twice, " \
                                 "you can ignore this warning."
            assert duplicated_warning in logged.output

    def test_fail_barcodes(self):
        barcode_fail_sheet = pathlib.Path(__file__).parent /"data" /"utilities_test" / "runsheet-barcode-fail.xlsx"
        sheet_data = pd.read_excel(barcode_fail_sheet, usecols = "A:C", skiprows = 3,
                                   dtype = {"KMA nr": str})
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nBarcodes ['RB3'] are not valid barcodes. " \
                    "Barcodes must consist of RB + a number between 01 and 96."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data, check_barcodes = True)

    def test_fail_no_barcodes(self):
        no_barcode_sheet = pathlib.Path(__file__).parent / "data" /"utilities_test" / "runsheet-no-barcode.xlsx"
        sheet_data = pd.read_excel(no_barcode_sheet, usecols = "A:C", skiprows = 3,
                                   dtype = {"KMA nr": str})
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nNo barcodes found."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data, check_barcodes = True)

    def test_fail_more_barcodes(self):
        more_barcodes_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-more-barcodes.xlsx"
        sheet_data = pd.read_excel(more_barcodes_sheet, usecols = "A:C", skiprows = 3,
                                   dtype = {"KMA nr": str})
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nAmount of sample IDs and barcodes don't match. " \
                    "There are 2 sample IDs but 3 barcodes."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data, check_barcodes = True)

    def test_fail_duplicated_barcodes(self):
        duplicated_barcodes_sheet = pathlib.Path(
            __file__).parent / "data" / "utilities_test" / "runsheet-barcode-duplication.xlsx"
        sheet_data = pd.read_excel(duplicated_barcodes_sheet, usecols = "A:C", skiprows = 3,
                                   dtype = {"KMA nr": str})
        error_msg = "The following issue(s) were detected with the runsheet:" \
                    "\nBarcode(s) ['RB02'] are duplicated."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            check_runsheet.check_sheet_format(sheet_data, check_barcodes = True)