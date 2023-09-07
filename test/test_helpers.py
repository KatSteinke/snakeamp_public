import pathlib
import re
import unittest

import pytest

import helpers


class TestFindRundir(unittest.TestCase):
    def test_find_absolute_path_success(self):
        test_path = pathlib.Path(__file__).parent / "data" / "helpers" \
                    / "test_dir_multi_pass" / "rawdata" / "subdir" / "fastq_pass"
        true_path = pathlib.Path(__file__).parent / "data" / "helpers" \
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
        with pytest.raises(FileNotFoundError, match = re.escape(error_msg)), \
                self.assertLogs("helpers", level = "INFO") as logged:
            log_msg = f"INFO:helpers:No barcode directories found in {test_path}.\n" \
                      f"Searching for barcodes in {test_path}/rawdata/*/fastq_pass..."
            helpers.get_fastq_pass_dir(test_path)
        assert log_msg in logged.output

    def test_fastq_fail_path(self):
        test_path = pathlib.Path(__file__).parent / "data" / "helpers" / "subdir_3" / "fastq_fail"
        error_msg = "fastq_pass folder(s) not found in expected location:\n" \
                    f"{test_path}/rawdata/*/fastq_pass\n" \
                    "Ensure correct directory and/or directory structure is used.\n" \
                    "Aborting 16S pipeline..."
        with pytest.raises(FileNotFoundError, match = re.escape(error_msg)),\
                self.assertLogs("helpers", level = "INFO") as logged:
            log_msg = f"INFO:helpers:No barcode directories found in {test_path}.\n" \
                      f"Searching for barcodes in {test_path}/rawdata/*/fastq_pass..."
            helpers.get_fastq_pass_dir(test_path)
        assert log_msg in logged.output

    def test_fail_multiple_fastq_dirs(self):
        test_path = pathlib.Path(__file__).parent / "data" / "helpers" / "test_dir_multi_pass"
        error_msg = f"The directory {test_path} contains " \
                    f"multiple fastq_pass directories." \
                    "Please choose the one containing the fastq files you want to analyze and " \
                    "specify the entire path to the fastq_pass directory."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            helpers.get_fastq_pass_dir(test_path)

    def test_fail_no_barcode_dirs(self):
        test_path = pathlib.Path(__file__).parent / "data" / "helpers" / \
                    "test_dir_multi_pass" / "rawdata" / "subdir_3" / "fastq_pass"
        error_msg = "fastq_pass or a subdirectory has been given " \
                    "but no barcode directories were found. " \
                    "Please give the path to the base directory " \
                    "or a directory containing barcode directories ('barcodeXX')."
        with pytest.raises(FileNotFoundError, match = re.escape(error_msg)):
            helpers.get_fastq_pass_dir(test_path)


class TestTranslateSampleNumbers(unittest.TestCase):
    number_to_letter = {"40": "H", "70": "P"}

    def test_wrong_sample_in(self):
        with pytest.raises(ValueError,
                           match="Invalid initial sample format int. Sample format can only be "
                                 "number or letter"):
            helpers.get_number_letter_combination(self.number_to_letter, "int", "letter")

    def test_wrong_sample_out(self):
        with pytest.raises(ValueError,
                           match="Invalid desired sample format str. Sample format can only be "
                                 "number or letter"):
            helpers.get_number_letter_combination(self.number_to_letter, "number", "str")

    def test_correct_results(self):
        number_to_number = helpers.get_number_letter_combination(self.number_to_letter, "number",
                                                                 "number")
        number_to_letter = helpers.get_number_letter_combination(self.number_to_letter, "number",
                                                                 "letter")
        letter_to_letter = helpers.get_number_letter_combination(self.number_to_letter, "letter",
                                                                 "letter")
        letter_to_number = helpers.get_number_letter_combination(self.number_to_letter, "letter",
                                                                 "number")
        self.assertEqual(number_to_number, {"40": "40", "70": "70"})
        self.assertEqual(number_to_letter, self.number_to_letter)
        self.assertEqual(letter_to_letter, {"H": "H", "P": "P"})
        self.assertEqual(letter_to_number, {"H": "40", "P": "70"})
