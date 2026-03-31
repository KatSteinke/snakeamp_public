import logging
import pathlib
import re
import unittest

import pandas as pd
import pytest

import helpers


class TestGetControls(unittest.TestCase):
    def test_no_controls(self):
        """Return unmatchable patterns for no controls."""
        expected_negative = re.compile('(?!.*)')
        expected_positive = re.compile('(?!.*)')
        test_negative, test_positive = helpers.get_control_patterns()
        assert expected_negative.pattern == test_negative.pattern
        assert expected_positive.pattern == test_positive.pattern

    def test_get_negative_control(self):
        """Return the correct pattern if there is a negative control."""
        expected_negative = re.compile('NegK')
        expected_positive = re.compile('(?!.*)')
        test_negative, test_positive = helpers.get_control_patterns(negative_control = "NegK")
        assert expected_negative.pattern == test_negative.pattern
        assert expected_positive.pattern == test_positive.pattern

    def test_get_positive_control(self):
        """Return the correct pattern if there is a positive control."""
        expected_negative = re.compile('(?!.*)')
        expected_positive = re.compile('PosK|PosK2')
        test_negative, test_positive = helpers.get_control_patterns(positive_control = {"PosK":
                                                                                            "Placeholderia bielefeldensis",
                                                                                        "PosK2":
                                                                                            "Placeholderia fakeorum"})
        assert expected_negative.pattern == test_negative.pattern
        assert expected_positive.pattern == test_positive.pattern

    def test_get_both_controls(self):
        """Return the correct pattern if both positive and negative controls are present."""
        expected_negative = re.compile('NegK')
        expected_positive = re.compile('PosK|PosK2')
        test_negative, test_positive = helpers.get_control_patterns(negative_control = "NegK",
                                                                    positive_control = {"PosK":
                                                                                            "Placeholderia bielefeldensis",
                                                                                        "PosK2":
                                                                                            "Placeholderia fakeorum"})
        assert expected_negative.pattern == test_negative.pattern
        assert expected_positive.pattern == test_positive.pattern


class TestGetPatterns(unittest.TestCase):
    def test_samples_only(self):
        """Create a sample number pattern without controls."""
        expected_pattern = re.compile(r"(test(?P<suffix>A|BC))")
        sample_number_pattern = "test(?P<suffix>A|BC)"
        test_pattern = helpers.get_id_pattern(sample_number_pattern)
        assert test_pattern == expected_pattern

    def test_fail_sample_number_blank(self):
        """Fail when there is no sample number pattern"""
        sample_number_pattern = ""
        error_msg = "Sample number pattern cannot be blank."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            helpers.get_id_pattern(sample_number_pattern)

    def test_negative_control(self):
        """Create a sample number pattern with a negative control."""
        expected_pattern = re.compile(r"(test(?P<suffix>A|BC)|NegK)")
        sample_number_pattern = "test(?P<suffix>A|BC)"
        negk_pattern = "NegK"
        test_pattern = helpers.get_id_pattern(sample_number_pattern,
                                              negative_control = negk_pattern)
        assert test_pattern == expected_pattern

    def test_positive_control(self):
        """Create a sample number pattern with a single positive control."""
        expected_pattern = re.compile(r"(test(?P<suffix>A|BC)|(PosK))")
        sample_number_pattern = "test(?P<suffix>A|BC)"
        positive_controls = {"PosK": "Placeholderia"}
        test_pattern = helpers.get_id_pattern(sample_number_pattern,
                                              positive_control = positive_controls)
        assert test_pattern == expected_pattern

    def test_multiple_positive_controls(self):
        """Create a sample number pattern with multiple positive controls."""
        expected_pattern = re.compile(r"(test(?P<suffix>A|BC)|(PosK|PosK2))")
        sample_number_pattern = "test(?P<suffix>A|BC)"
        positive_controls = {"PosK": "Placeholderia", "PosK2": "Fakeobacter"}
        test_pattern = helpers.get_id_pattern(sample_number_pattern,
                                              positive_control = positive_controls)
        assert test_pattern == expected_pattern

    def test_all_control_types(self):
        """Create a sample number pattern with both positive and negative controls."""
        expected_pattern = re.compile(r"(test(?P<suffix>A|BC)|NegK|(PosK|PosK2))")
        sample_number_pattern = "test(?P<suffix>A|BC)"
        negative_control = "NegK"
        positive_controls = {"PosK": "Placeholderia", "PosK2": "Fakeobacter"}
        test_pattern = helpers.get_id_pattern(sample_number_pattern,
                                              negative_control = negative_control,
                                              positive_control = positive_controls)
        assert test_pattern == expected_pattern


class TestCheckBarcodeDirs(unittest.TestCase):
    def test_raise_when_no_barcodes(self):
        """Raise an error if no barcode directories are present."""
        test_path = pathlib.Path(__file__).parent / "data" / "helpers" / \
                    "test_dir_multi_pass" / "rawdata" / "subdir_3" / "fastq_pass"
        error_msg = "No barcode directories found in fastq_pass directory."
        with pytest.raises(FileNotFoundError, match = re.escape(error_msg)):
            helpers.check_barcode_dirs(test_path)

class TestFindRundir(unittest.TestCase):
    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    def test_find_absolute_path_success(self):
        test_path = pathlib.Path(__file__).parent / "data" / "helpers" \
                    / "test_dir_multi_pass" / "rawdata" / "subdir" / "fastq_pass"
        true_path = pathlib.Path(__file__).parent / "data" / "helpers" \
                    / "test_dir_multi_pass" / "rawdata" / "subdir"
        log_msg = "Data is retrieved from the following folder:\n" \
                  f"{true_path}"
        with self._caplog.at_level(logging.INFO, logger = "helpers"):
            test_fastq = helpers.get_fastq_pass_parent(test_path)
            assert ("helpers", logging.INFO, log_msg) in self._caplog.record_tuples
        assert test_fastq == true_path

    def test_handle_parent_dir(self):
        test_path = pathlib.Path(__file__).parent / "data" / "helpers" \
                    / "test_dir_multi_pass" / "rawdata" / "subdir"
        true_path = pathlib.Path(__file__).parent / "data" / "helpers" \
                    / "test_dir_multi_pass" / "rawdata" / "subdir"
        log_msg = "Data is retrieved from the following folder:\n" \
                  f"{true_path}"
        with self._caplog.at_level(logging.INFO, logger = "helpers"):
            test_fastq = helpers.get_fastq_pass_parent(test_path)
            assert ("helpers", logging.INFO, log_msg) in self._caplog.record_tuples
        assert test_fastq == true_path

    def test_return_fastq_when_given(self):
        test_path = pathlib.Path(__file__).parent / "data" / "snake_helpers" / "test_dir"
        true_path = pathlib.Path(__file__).parent / "data" / "snake_helpers" / "test_dir" \
                    / "rawdata" / "subdir"
        log_msg = "Data is retrieved from the following folder:\n" \
                  f"{true_path}"
        with self._caplog.at_level(logging.INFO, logger = "helpers"):
            test_fastq = helpers.get_fastq_pass_parent(test_path)
            assert ("helpers", logging.INFO, log_msg) in self._caplog.record_tuples
        assert test_fastq == true_path

    def test_handle_different_samples(self):
        """Handle different sample names."""
        test_path = pathlib.Path(__file__).parent / "data" / "helpers" / "test_dir"
        true_path = pathlib.Path(__file__).parent / "data" / "helpers" / "test_dir" \
                    / "no_sample" / "subdir"
        log_msg = "Data is retrieved from the following folder:\n" \
                  f"{true_path}"
        with self._caplog.at_level(logging.INFO, logger = "helpers"):
            test_fastq = helpers.get_fastq_pass_parent(test_path)
            assert ("helpers", logging.INFO, log_msg) in self._caplog.record_tuples
        assert test_fastq == true_path

    def test_fail_path(self):
        test_path = pathlib.Path(__file__).parent / "data" / "helpers" / "subdir"
        error_msg = f"fastq_pass folder not found in {test_path} or any subfolders. \n" \
                    "Ensure correct directory and/or directory structure is used.\n" \
                    "Aborting pipeline..."
        log_msg = f"Searching for fastq_pass folder in {test_path}..."
        with pytest.raises(FileNotFoundError, match = re.escape(error_msg)), \
                self._caplog.at_level(logging.INFO, logger = "helpers"):
            helpers.get_fastq_pass_parent(test_path)
        assert ("helpers", logging.INFO, log_msg) in self._caplog.record_tuples

    def test_fastq_fail_path(self):
        test_path = pathlib.Path(__file__).parent / "data" / "helpers" / "subdir_3" / "fastq_fail"
        error_msg = f"fastq_pass folder not found in {test_path} or any subfolders. \n" \
                    "Ensure correct directory and/or directory structure is used.\n" \
                    "Aborting pipeline..."
        log_msg = f"Searching for fastq_pass folder in {test_path}..."
        with pytest.raises(FileNotFoundError, match = re.escape(error_msg)), \
                self._caplog.at_level(logging.INFO, logger = "helpers"):
            helpers.get_fastq_pass_parent(test_path)
        assert ("helpers", logging.INFO, log_msg) in self._caplog.record_tuples

    def test_fail_multiple_fastq_dirs(self):
        test_path = pathlib.Path(__file__).parent / "data" / "helpers" / "test_dir_multi_pass"
        error_msg = f"The directory {test_path} contains " \
                    f"multiple fastq_pass directories." \
                    "Please choose the one containing the fastq files you want to analyze and " \
                    "specify the entire path to the fastq_pass directory."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            helpers.get_fastq_pass_parent(test_path)

    def test_warn_no_barcode_dirs(self):
        test_path = pathlib.Path(__file__).parent / "data" / "helpers" / \
                    "test_dir_multi_pass" / "rawdata" / "subdir_3" / "fastq_pass"
        warn_msg = "No barcode directories found in fastq_pass directory. " \
                   "This may be due to a delay in copying files from the sequencer, but could" \
                   " also mean you have given the wrong path. \n" \
                   "Only continue if you are sure. "
        with self._caplog.at_level(logging.WARNING, logger = "helpers"):
            helpers.get_fastq_pass_parent(test_path)
            assert ("helpers", logging.WARNING, warn_msg) in self._caplog.record_tuples

    def test_fail_nonexistent_dir(self):
        """Fail if a nonexistent directory is supplied."""
        good_path = (pathlib.Path(__file__).parent / "data" / "helpers" / "scenario3"
                     / "rawdata" / "subdir")
        test_find = helpers.get_fastq_pass_parent(good_path)
        assert good_path == test_find
        test_path = (pathlib.Path(__file__).parent / "data" / "helpers" / "scenario3"
                     / "subdir" / "fastq_pass")
        error_msg =  f"The supplied folder {test_path} does not exist. \n" \
                    "Ensure correct directory and/or directory structure is used.\n" \
                    "Aborting 16S pipeline..."
        with pytest.raises(FileNotFoundError, match = re.escape(error_msg)):
            helpers.get_fastq_pass_parent(test_path)


class TestTranslateSampleType(unittest.TestCase):
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
        number_format = re.compile(r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})')
        positive_control = re.compile('PosK')
        negative_control = re.compile('NegK')
        expected_prefix = "11"
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
        number_format = re.compile(r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})')
        positive_control = re.compile('PosK')
        negative_control = re.compile('NegK')
        test_prefix = helpers.extract_sample_number_part(test_number, "sample_type",
                                                         number_format,
                                                         negative_control, positive_control)
        assert pd.isna(test_prefix)

    def test_handle_negative_control(self):
        """Don't try to extract the component from a negative control."""
        test_number = "NegK"
        number_format = re.compile(r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})')
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

class TestRearrangeSampleNumber(unittest.TestCase):
    def test_no_change(self):
        input_number = "F99123456-1"
        true_number = "F99123456-1"
        pattern_in = re.compile(r"(?P<sample_type>[BDFT])"
                                r"(?P<sample_year>\d{2})"
                                r"(?P<sample_number>\d{6})"
                                r"(?P<bact_number>-\d)")
        order_out = {1:"sample_type", 2: "sample_year", 3: "sample_number",4: "bact_number"}
        test_number = helpers.rearrange_sample_number(input_number, pattern_in, order_out)
        assert test_number == true_number

    def test_rearrange_success(self):
        input_number = "F99123456-1"
        true_number = "F12345699-1"
        pattern_in = re.compile(r"(?P<sample_type>[BDFT])"
                                r"(?P<sample_year>\d{2})"
                                r"(?P<sample_number>\d{6})"
                                r"(?P<bact_number>-\d)")
        order_out ={1:"sample_type", 2: "sample_number", 3: "sample_year",4: "bact_number"}
        test_number = helpers.rearrange_sample_number(input_number, pattern_in, order_out)
        assert test_number == true_number

    def test_fail_too_many_components(self):
        input_number = "F99123456-1"
        pattern_in = re.compile(r"(?P<sample_type>[BDFT])"
                                r"(?P<sample_number>\d{6})"
                                r"(?P<bact_number>-\d)")
        order_out = {1:"sample_type", 2: "sample_year", 3: "sample_number",4: "bact_number"}
        error_msg = "Not all desired sample number components could be " \
                    "found in the original format. " \
                    "Desired sample number format contains additional " \
                    "components {'sample_year'}."
        with pytest.raises(KeyError, match = error_msg):
            helpers.rearrange_sample_number(input_number, pattern_in, order_out)

    def test_fail_no_hits(self):
        input_number = "1199123456-1"
        pattern_in = re.compile(r"(?P<sample_type>[BDFT])"
                                r"(?P<sample_year>\d{2})"
                                r"(?P<sample_number>\d{6})"
                                r"(?P<bact_number>-\d)")
        order_out = {1:"sample_type", 2: "sample_year", 3: "sample_number",4: "bact_number"}
        error_msg = "No match in sample number 1199123456-1."
        with pytest.raises(ValueError, match = error_msg):
            helpers.rearrange_sample_number(input_number, pattern_in, order_out)

    def test_success_remove_component(self):
        input_number = "F99123456-1"
        true_number = "F99123456"
        pattern_in = re.compile(r"(?P<sample_type>[BDFT])"
                                r"(?P<sample_year>\d{2})"
                                r"(?P<sample_number>\d{6})"
                                r"(?P<bact_number>-\d)")
        order_out = {1:"sample_type", 2: "sample_year", 3: "sample_number"}
        test_number = helpers.rearrange_sample_number(input_number, pattern_in, order_out)
        assert test_number == true_number


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


class TestTranslateSampleNumber(unittest.TestCase):
    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    def test_fail_no_match(self):
        """Complain if the input sample number does not match the original format."""
        format_in = re.compile('(?P<sample_type>[BDPT]|[1357]0)(?P<sample_number>\d{8})(?P<bact_number>-\d)')
        format_out = re.compile('(?P<sample_type>[BDPT]|[1357]0)(?P<sample_number>\d{8})(?P<bact_number>-\d)')
        prefix_mapping = { "70": "P", "30": "B", "10": "D", "50": "T"}
        sample_number = "F99123456-1"
        error_msg = "Sample number F99123456-1 does not match specified input format."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            helpers.translate_sample_number(sample_number, format_in, format_out, prefix_mapping,
                                            re.compile("PosK"), re.compile('NegK[a-zA-Z0-9_-]*'))

    def test_fail_missing_prefix(self):
        """Complain if the sample number's prefix is not contained in the prefix mapping."""
        format_in = re.compile(
            '(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_number>\d{8})(?P<bact_number>-\d)')
        format_out = re.compile(
            '(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_number>\d{8})(?P<bact_number>-\d)')
        prefix_mapping = {"70": "P", "30": "B", "10": "D", "11": "F", "50": "T"}
        sample_number = "F99123456-1"
        error_msg = "Prefix F not found (allowed prefixes are ['70', '30', '10', '11', '50'])."
        with pytest.raises(KeyError, match = re.escape(error_msg)):
            helpers.translate_sample_number(sample_number, format_in, format_out, prefix_mapping,
                                            re.compile("PosK"), re.compile('NegK[a-zA-Z0-9_-]*'))

    def test_fail_no_match_after_translate(self):
        """Complain if the sample number does not match the desired format after translation."""
        format_in = re.compile(
            '(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_number>\d{8})(?P<bact_number>-\d)')
        format_out = re.compile(
            '(?P<sample_type>[1357]0|11)(?P<sample_number>\d{8})(?P<bact_number>-\d)')
        prefix_mapping = {"70": "P", "30": "B", "10": "D", "11": "F", "50": "T"}
        sample_number = "1199123456-1"
        error_msg = ("Translated sample number F99123456-1 (was 1199123456-1)"
                     " does not match desired output format.")
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            helpers.translate_sample_number(sample_number, format_in, format_out, prefix_mapping,
                                            re.compile("PosK"), re.compile('NegK[a-zA-Z0-9_-]*'))

    def test_translate_no_change(self):
        """Pass the sample number through without any changes if none are needed."""
        format_in = re.compile(
            '(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_number>\d{8})(?P<bact_number>-\d)')
        format_out = re.compile(
            '(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_number>\d{8})(?P<bact_number>-\d)')
        prefix_mapping = {"P": "P", "B": "B", "D": "D", "F": "F", "T": "T"}
        sample_number = "F99123456-1"
        test_number = helpers.translate_sample_number(sample_number, format_in, format_out,
                                                      prefix_mapping, re.compile("PosK"),
                                                      re.compile('NegK[a-zA-Z0-9_-]*'))
        assert sample_number == test_number

    def test_translate_prefix_only(self):
        """Translate a prefix."""
        format_in = re.compile(
            '(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_number>\d{8})(?P<bact_number>-\d)')
        format_out = re.compile(
            '(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_number>\d{8})(?P<bact_number>-\d)')
        prefix_mapping = {"70": "P", "30": "B", "10": "D", "11": "F", "50": "T"}
        sample_number = "1199123456-1"
        expected_number = "F99123456-1"
        test_number = helpers.translate_sample_number(sample_number, format_in, format_out,
                                                      prefix_mapping, re.compile("PosK"),
                                                      re.compile('NegK[a-zA-Z0-9_-]*'))
        assert expected_number == test_number

    def test_rearrange_only(self):
        """Rearrange a sample number."""
        format_in = re.compile(
            '(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_number>\d{8})(?P<bact_number>-\d)')
        format_out = re.compile(
            '(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_number>\d{8})')
        prefix_mapping = {"P": "P", "B": "B", "D": "D", "F": "F", "T": "T"}
        sample_number = "F99123456-1"
        expected_number = "F99123456"
        test_number = helpers.translate_sample_number(sample_number, format_in, format_out,
                                                      prefix_mapping, re.compile("PosK"),
                                                      re.compile('NegK[a-zA-Z0-9_-]*'))
        assert expected_number == test_number

    def test_translate_and_rearrange(self):
        """Translate the prefix and rearrange the sample number."""
        format_in = re.compile(
            '(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_number>\d{8})(?P<bact_number>-\d)')
        format_out = re.compile(
            '(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_number>\d{8})')
        prefix_mapping = {"70": "P", "30": "B", "10": "D", "11": "F", "50": "T"}
        sample_number = "1199123456-1"
        expected_number = "F99123456"
        test_number = helpers.translate_sample_number(sample_number, format_in, format_out,
                                                      prefix_mapping, re.compile("PosK"),
                                                      re.compile('NegK[a-zA-Z0-9_-]*'))
        assert expected_number == test_number

    def test_handle_controls(self):
        """Don't translate positive or negative controls"""
        format_in = re.compile(
            '(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_number>\d{8})(?P<bact_number>-\d)')
        format_out = re.compile(
            '(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_number>\d{8})')
        prefix_mapping = {"70": "P", "30": "B", "10": "D", "11": "F", "50": "T"}
        negative_control = "NegK"
        test_negative_control = helpers.translate_sample_number(negative_control, format_in,
                                                                format_out,
                                                                prefix_mapping, re.compile("PosK"),
                                                                re.compile('NegK[a-zA-Z0-9_-]*'))
        assert negative_control == test_negative_control
        positive_control = "PosK"
        test_positive_control = helpers.translate_sample_number(positive_control, format_in,
                                                                format_out, prefix_mapping,
                                                                re.compile("PosK"),
                                                                re.compile('NegK[a-zA-Z0-9_-]*'))
        assert positive_control == test_positive_control

    def test_handle_no_type(self):
        """Handle a sample number format without sample type."""
        format_in = re.compile(
            '(?P<sample_number>\d{8})(?P<bact_number>-\d)')
        format_out = re.compile(
            '(?P<sample_number>\d{8})')
        prefix_mapping = {"70": "P", "30": "B", "10": "D", "11": "F", "50": "T"}
        sample_number = "99123456-1"
        expected_number = "99123456"
        log_msg = ("No sample_type given in sample number format specification;"
                   " cannot translate sample type.")

        with self._caplog.at_level(logging.INFO, logger = "helpers"):
            test_number = helpers.translate_sample_number(sample_number, format_in, format_out,
                                                          prefix_mapping, re.compile("PosK"),
                                                          re.compile('NegK[a-zA-Z0-9_-]*'))
            assert ("helpers", logging.INFO, log_msg) in self._caplog.record_tuples
        assert expected_number == test_number


class TestCheckExperimentName(unittest.TestCase):
    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    def test_no_issues(self):
        """Ensure a name without problematic components doesn't raise an exception."""
        test_name = "samplerun"
        log_msg = "No issues found with experiment name samplerun."
        with self._caplog.at_level(logging.DEBUG, logger = "helpers"):
            helpers.check_experiment_name_problems(test_name)
        assert ("helpers", logging.DEBUG, log_msg) in self._caplog.record_tuples

    def test_check_whitespace_breaks(self):
        """Ensure that the script complains on names containing whitespace"""
        test_name = "sample run"
        with pytest.raises(ValueError, match= "The run name contains spaces or line breaks."):
            helpers.check_experiment_name_problems(test_name)

    def test_check_illegal_in_windows(self):
        """Ensure that the script complains for names that are illegal in Windows"""
        test_control_char = "NUL"
        test_bad_char = "sample:run"
        with pytest.raises(ValueError,
                           match = "The run name contains a character that cannot be used "
                                   "in Windows filenames."):
            helpers.check_experiment_name_problems(test_control_char)
        with pytest.raises(ValueError,
                           match = "The run name contains a character that cannot be used"
                                   " in Windows filenames."):
            helpers.check_experiment_name_problems(test_bad_char)

    def test_check_slash_breaks(self):
        """Ensure the script complains for names containing forward slashes"""
        test_name = "sample/run"
        error_msg = "The run name contains a forward slash (/). " \
                    "This will break the result directory." \
                    " Replace forward slashes with underscores (_)."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            helpers.check_experiment_name_problems(test_name)


class TestExtractNanoporeRun(unittest.TestCase):
    def test_get_run_name(self):
        test_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                     / "test_nanopore_runsheet.xlsx"
        true_run_name = "NANO_Amplicon_Y20990101_RUN0001_XYZ"
        test_run_name = helpers.extract_nanopore_run_name(test_sheet)
        assert test_run_name == true_run_name

    def test_get_short_run_name(self):
        """Extract a run name without zero padding."""
        test_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                     / "test_nanopore_runsheet_short_name.xlsx"
        true_run_name = "NANO_Amplicon_Y20990101_RUN1_XYZ"
        test_run_name = helpers.extract_nanopore_run_name(test_sheet)
        assert test_run_name == true_run_name

    def test_check_whitespace_breaks(self):
        test_sheet = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                     / "nanopore_bad_name.xlsx"
        with pytest.raises(ValueError, match= "The run name contains spaces or line breaks."):
            helpers.extract_nanopore_run_name(test_sheet)

    def test_check_illegal_in_windows(self):
        test_control_char = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                     / "nanopore_controlchar.xlsx"
        test_bad_char = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                     / "nanopore_badchar.xlsx"
        with pytest.raises(ValueError,
                           match = "The run name contains a character that cannot be used "
                                   "in Windows filenames."):
            helpers.extract_nanopore_run_name(test_control_char)
        with pytest.raises(ValueError,
                           match = "The run name contains a character that cannot be used"
                                   " in Windows filenames."):
            helpers.extract_nanopore_run_name(test_bad_char)

    def test_slash_breaks(self):
        test_with_slash = pathlib.Path(__file__).parent / "data" / "utilities_test" \
                     / "nanopore_slash.xlsx"
        error_msg = "The run name contains a forward slash (/). " \
                    "This will break the result directory." \
                    " Replace forward slashes with underscores (_)."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            helpers.extract_nanopore_run_name(test_with_slash)