import pathlib
import re
import unittest

import pandas as pd
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


class TestFindPart(unittest.TestCase):
    def test_get_match(self):
        """Ensure a component matching the pattern is reported."""
        test_number = "1199123456"
        number_format = re.compile(r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})')
        positive_control = re.compile('PosK')
        negative_control = re.compile('NegK')
        expected_prefix = "70"
        test_prefix = helpers.extract_sample_number_part(test_number, "sample_type",
                                                         number_format,
                                                         negative_control, positive_control)
        assert test_prefix == expected_prefix
        expected_year = "99"
        test_year = helpers.extract_sample_number_part(test_number, "sample_year",
                                                       number_format,
                                                       negative_control, positive_control)
        assert test_year == expected_year

    def test_handle_positive_control(self):
        """Don't try to extract the component from a positive control."""
        test_number = "PosK"
        number_format = re.compile(r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})')
        positive_control = re.compile('PosK')
        negative_control = re.compile('NegK')
        test_prefix = helpers.extract_sample_number_part(test_number, "sample_type",
                                                         number_format,
                                                         negative_control, positive_control)
        assert pd.isna(test_prefix)

    def test_handle_negative_control(self):
        """Don't try to extract the component from a negative control."""
        test_number = "NegK"
        number_format = re.compile(r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})')
        positive_control = re.compile('PosK')
        negative_control = re.compile('NegK')
        test_prefix = helpers.extract_sample_number_part(test_number, "sample_type",
                                                         number_format,
                                                         negative_control, positive_control)
        assert pd.isna(test_prefix)

    def test_handle_no_hit(self):
        """Handle a number not matching the pattern."""
        test_number = "11123456"
        number_format = re.compile(r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})')
        positive_control = re.compile('PosK')
        negative_control = re.compile('NegK')
        test_prefix = helpers.extract_sample_number_part(test_number, "sample_type",
                                                         number_format,
                                                         negative_control, positive_control)
        assert pd.isna(test_prefix)

    def test_warn_no_prefix(self):
        """Alert the user if the pattern doesn't contain a definition of the component."""
        test_number = "99123456"
        number_format = re.compile(r'(?P<sample_year>\d{2})(?P<sample_number>\d{6})')
        positive_control = re.compile('PosK')
        negative_control = re.compile('NegK')
        log_msg = "WARNING:helpers:Group name sample_type not found in sample number pattern." \
                  " Component cannot be extracted."
        with self.assertLogs("helpers") as logged:
            test_prefix = helpers.extract_sample_number_part(test_number, "sample_type",
                                                             number_format,
                                                             negative_control, positive_control)
            assert log_msg in logged.output
        assert pd.isna(test_prefix)


class TestExtractMatchGroup(unittest.TestCase):
    def test_extract_simple_group(self):
        """Extract a simple match group."""
        sample_pattern = re.compile(r"something(?P<testgroup>\d)somethingelse")
        expected_pattern = re.compile(r"(?P<testgroup>\d)")
        test_pattern = helpers.parse_out_group_pattern(sample_pattern, "testgroup")
        assert test_pattern == expected_pattern

    def test_extract_group_internal_parentheses(self):
        """Extract a group with internal parentheses."""
        sample_pattern = re.compile(r"something(?P<testgroup>\d(AB|CD))somethingelse")
        expected_pattern = re.compile(r"(?P<testgroup>\d(AB|CD))")
        test_pattern = helpers.parse_out_group_pattern(sample_pattern, "testgroup")
        assert test_pattern == expected_pattern

    def test_extract_one_only(self):
        """Extract one group from a pattern with multiple groups."""
        sample_pattern = re.compile(r"something(?P<testgroup>\d)(?P<test2>somethingelse)")
        expected_pattern = re.compile(r"(?P<testgroup>\d)")
        test_pattern = helpers.parse_out_group_pattern(sample_pattern, "testgroup")
        assert test_pattern == expected_pattern

    def test_extract_outer_nested(self):
        """Extract the outer group of a nested named group."""
        sample_pattern = re.compile(r"something(?P<testgroup>\d(?P<test2>XYZ))somethingelse")
        expected_pattern = re.compile(r"(?P<testgroup>\d(?P<test2>XYZ))")
        test_pattern = helpers.parse_out_group_pattern(sample_pattern, "testgroup")
        assert test_pattern == expected_pattern

    def test_extract_inner_nested(self):
        """Extract the inner group of a nested named group."""
        sample_pattern = re.compile(r"something(?P<testgroup>\d(?P<test2>XYZ))somethingelse")
        expected_pattern = re.compile(r"(?P<test2>XYZ)")
        test_pattern = helpers.parse_out_group_pattern(sample_pattern, "test2")
        assert test_pattern == expected_pattern

    def test_fail_group_not_found(self):
        """Complain if the pattern does not contain the group."""
        sample_pattern = re.compile(r"something(?P<testgroup>\d)somethingelse")
        error_msg = "Group test_group not found in named groups."
        with pytest.raises(KeyError, match = re.escape(error_msg)):
            helpers.parse_out_group_pattern(sample_pattern, "test_group")


class TestAddYearsInSheet(unittest.TestCase):
    def test_success_add_year(self):
        """Ensure year is added to properly formatted sample numbers."""
        test_input = pd.DataFrame(data = {"KMA nr": ["11123456", "11123456"],
                                          "årstal": ["99", "99"],
                                          "Barkode NB": ["RB01", "RB02"]})
        expected_df = pd.DataFrame(data = {"KMA nr": ["11123456", "11123456"],
                                           "årstal": ["99", "99"],
                                           "Barkode NB": ["RB01", "RB02"],
                                           "prøvenr": ["1199123456", "1199123456"]})
        test_df = helpers.add_years_in_sheet(test_input)
        pd.testing.assert_frame_equal(expected_df, test_df)

    def test_complain_no_year_col(self):
        """Ensure an error is raised if there is no column for the sample year."""
        test_input = pd.DataFrame(data = {"KMA nr": ["11123456", "11123456"],
                                          "Barkode NB": ["RB01", "RB02"]})
        error_msg = "No year column found. " \
                    "The pipeline needs a column named 'årstal' to add year to sample number."
        with pytest.raises(KeyError, match=re.escape(error_msg)):
            helpers.add_years_in_sheet(test_input)

    def test_complain_wrong_year_format(self):
        """Ensure an error is raised if the year is given in the wrong format."""
        test_input = pd.DataFrame(data = {"KMA nr": ["11123456", "11123456"],
                                          "årstal": ["99", "2099"],
                                          "Barkode NB": ["RB01", "RB02"]})
        bad_years = pd.DataFrame(data={"KMA nr": ["11123456"],
                                       "årstal": ["2099"],
                                       "Barkode NB": ["RB02"]})
        error_msg = "Invalid year values detected. Year must be given as YY only." \
                    " Affected samples:\n" \
                    f"{bad_years.to_string()}"
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            helpers.add_years_in_sheet(test_input)

    def test_complain_blank_year_column(self):
        """Ensure an error is raised if no year is given for a sample."""
        test_input = pd.DataFrame(data = {"KMA nr": ["11123456", "11123456"],
                                          "årstal": ["99", pd.NA],
                                          "Barkode NB": ["RB01", "RB02"]})
        bad_years = pd.DataFrame(data = {"KMA nr": ["11123456"],
                                         "årstal": [pd.NA],
                                         "Barkode NB": ["RB02"]})
        error_msg = "No year given for one or more samples. Please add a year to these samples." \
                    " Affected samples:\n" \
                    f"{bad_years.to_string()}"
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            helpers.add_years_in_sheet(test_input)

    def test_handle_controls(self):
        """Ensure positive and negative controls are processed unaltered."""
        test_input = pd.DataFrame(data = {"KMA nr": ["11123456", "11123456", "PosK", "NegK"],
                                          "årstal": ["99", "99", "", ""],
                                          "Barkode NB": ["RB01", "RB02", "RB03", "RB04"]})
        expected_df = pd.DataFrame(data = {"KMA nr": ["11123456", "11123456", "PosK", "NegK"],
                                           "årstal": ["99", "99", "", ""],
                                           "Barkode NB": ["RB01", "RB02", "RB03", "RB04"],
                                           "prøvenr": ["1199123456", "1199123456", "PosK", "NegK"]})
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
                                                      {"splice_in_date": True,
                                                       "length_without_date": 8,
                                                       "splice_after": 2},
                                                  "negative_control": 'NegK',
                                                  "positive_control": {"PosK": "Placeholderia"}},
                       "barcode_format": "RB[0-9]{2}",  # format of barcodes in runsheet
                       "barcode_prefix": "RB"  # barcode prefix as letter (for transferring original fastqs by barcode)
                       }
        test_df = helpers.add_years_in_sheet(test_input, active_config = test_config)
        pd.testing.assert_frame_equal(expected_df, test_df)
