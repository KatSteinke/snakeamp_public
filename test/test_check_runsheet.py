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

    def test_prefix_not_in_sheet(self):
        with pytest.raises(ValueError, match="No samples with prefix X found in runsheet."):
            check_runsheet.check_by_prefix(self.sheet_data, self.lab_info_data, "X", "P")

    def test_sample_not_in_report(self):
        fail_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test"\
                        / "test_notinmads_runsheet.xlsx"
        fail_data = pd.read_excel(fail_runsheet, usecols="A:C", skiprows=3,
                                  dtype={"KMA nr": str, "Barkode NB": str})
        fail_data = fail_data.dropna()
        error_msg = "Samples ['11410000'] were not found in MADS report. " \
                    "Please check that sample numbers are correct."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_by_prefix(fail_data, self.lab_info_data, "30", "B")

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
            check_runsheet.check_by_prefix(self.sheet_data, lab_info_data, "30", "B")


class TestCheckRunsheet(unittest.TestCase):
    def test_multiple_fails(self):
        fail_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_notinmads_runsheet.xlsx"
        fake_mads = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "fake_mads_data.csv"
        sheet_data = pd.read_excel(fail_runsheet, usecols="A:B", skiprows=3,
                                   dtype={"KMA nr": str})
        sheet_data = sheet_data.dropna()
        error_msg = "The following issues were encountered:\n" \
                    "Samples ['11400000'] were not found in MADS report. " \
                    "Please check that sample numbers are correct.\n" \
                    "Samples ['F0410000'] were not found in MADS report. " \
                    "Please check that sample numbers are correct."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_runsheet(sheet_data, fake_mads)

    def test_success(self):
        test_runsheet = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "test_translate_runsheet.xlsx"
        fake_mads = pathlib.Path(__file__).parent / "data" / "sample_sheet_test" / "fake_mads_data.csv"
        sheet_data = pd.read_excel(test_runsheet, usecols = "A:B", skiprows = 3,
                                   dtype = {"KMA nr": str})
        sheet_data = sheet_data.dropna()
        success_msg = "INFO:check_runsheet:The runsheet is correct."
        with self.assertLogs("check_runsheet") as logged:
            check_runsheet.check_runsheet(sheet_data, fake_mads)
            assert success_msg in logged.output


class TestCheckSampleNumbers(unittest.TestCase):
    def test_fail_ids(self):
        id_fail_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-id-fail.xlsx"
        error_msg = "Sample IDs ['123'] are not valid. " \
                    "Sample IDs must start with 70 or 30 or 10 or 50 followed by eight numbers" \
                    " (six if leaving out year). " \
                    "Negative controls must be given in the format . " \
                    "Please correct sample IDs in runsheet."
        sheet_data = pd.read_excel(id_fail_sheet, usecols="A:B", skiprows=3,
                                   dtype={"KMA nr": str})
        sheet_data = sheet_data.dropna()
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_sample_numbers(sheet_data)

    def test_fail_no_ids(self):
        no_id_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-no-id.xlsx"
        sheet_data = pd.read_excel(no_id_sheet, usecols = "A:B", skiprows = 3,
                                   dtype = {"KMA nr": str})
        sheet_data = sheet_data.dropna()
        error_msg = "No sample IDs found."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            check_runsheet.check_sample_numbers(sheet_data)

    def test_fail_no_positive_control(self):
        no_positive_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-no-posk.xlsx"
        sheet_data = pd.read_excel(no_positive_sheet, usecols="A:B", skiprows=3,
                                   dtype={"KMA nr": str})
        sheet_data = sheet_data.dropna()
        error_msg = "No positive controls given in runsheet."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_sample_numbers(sheet_data)

    def test_fail_no_negative_control(self):
        no_negative_sheet = pathlib.Path(__file__).parent / "data" /"utilities_test" / "runsheet-no-negk.xlsx"
        sheet_data = pd.read_excel(no_negative_sheet, usecols="A:B", skiprows=3,
                                   dtype={"KMA nr": str})
        sheet_data = sheet_data.dropna()
        error_msg = "No negative controls given in runsheet."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_runsheet.check_sample_numbers(sheet_data)

    def test_warn_duplicated_ids(self):
        duplicated_id_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" / "runsheet-id-duplication.xlsx"
        sheet_data = pd.read_excel(duplicated_id_sheet, usecols="A:B", skiprows=3,
                                   dtype={"KMA nr": str})
        with self.assertLogs("check_runsheet") as logged:
            check_runsheet.check_sample_numbers(sheet_data)
            duplicated_warning = "WARNING:check_runsheet:Sample number(s) ['1123456789'] are duplicated." \
                                 " If you are sure you want to sequence the same sample twice, " \
                                 "you can ignore this warning."
            assert duplicated_warning in logged.output
