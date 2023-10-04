import pathlib
import re
import unittest

import numpy as np
import pandas as pd
import pytest

import summarize_emu


class TestExtractCounts(unittest.TestCase):
    workflow_config = {"sample_number_settings": {"sample_number_format": r"barcode\d{2}",
                                                  "positive_control": {},
                                                  "negative_control": ""},
                       "barcode_format": "RB[0-9]{2}"}

    def test_get_counts_success(self):
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "barcode01_RB01_rel-abundance.tsv"
        expected_results = pd.DataFrame(data={"abundance_from_all": [0.75, 0.2, 0.05],
                                              "estimated counts": [15.0, 4.0, 1.0],
                                              "medtages": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia fakeorum",
                                                                 "Placeholderia bielefeldensis",
                                                                 "unassigned"], name = "species"))
        barcode_header = ["barcode01_RB01"] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([barcode_header,
                                                              expected_results.columns])
        test_results = summarize_emu.report_species_per_barcode(sample_path, self.workflow_config)
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_handle_isolate_number(self):
        """Handle formats including special characters"""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      # TODO: replace with 0-9
                                                      "positive_control": {},
                                                      "negative_control": ""},
                           "barcode_format": "RB[0-9]{2}"}
        sample_path = pathlib.Path(__file__).parent / "data" / "summarize_emu" / "F99123456-0_RB01_rel-abundance.tsv"
        expected_results = pd.DataFrame(data = {"abundance_from_all": [0.75, 0.2, 0.05],
                                                "estimated counts": [15.0, 4.0, 1.0],
                                                "medtages": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia fakeorum",
                                                                 "Placeholderia bielefeldensis",
                                                                 "unassigned"], name = "species"))
        barcode_header = ["F99123456-0_RB01"] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([barcode_header,
                                                              expected_results.columns])
        test_results = summarize_emu.report_species_per_barcode(sample_path, workflow_config)
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_handle_negative_control(self):
        """Handle negative controls"""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      # TODO: replace with 0-9
                                                      "positive_control": {},
                                                      "negative_control": "NegK"},
                       "barcode_format": "RB[0-9]{2}"}
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "NegK_RB02_rel-abundance.tsv"
        expected_results = pd.DataFrame(data = {"abundance_from_all": [0.75, 0.2, 0.05],
                                                "estimated counts": [15.0, 4.0, 1.0],
                                                "medtages": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia fakeorum",
                                                                 "Placeholderia bielefeldensis",
                                                                 "unassigned"], name = "species"))
        barcode_header = ["NegK_RB02"] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([barcode_header,
                                                              expected_results.columns])
        test_results = summarize_emu.report_species_per_barcode(sample_path, workflow_config)
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_handle_positive_control(self):
        """Handle positive controls"""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      # TODO: replace with 0-9
                                                      "positive_control": {"PosK": "Placeholderia"},
                                                      "negative_control": ""},
                       "barcode_format": "RB[0-9]{2}"}
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "PosK_RB03_rel-abundance.tsv"
        expected_results = pd.DataFrame(data = {"abundance_from_all": [0.75, 0.2, 0.05],
                                                "estimated counts": [15.0, 4.0, 1.0],
                                                "medtages": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia fakeorum",
                                                                 "Placeholderia bielefeldensis",
                                                                 "unassigned"], name = "species"))
        barcode_header = ["PosK_RB03"] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([barcode_header,
                                                              expected_results.columns])
        test_results = summarize_emu.report_species_per_barcode(sample_path, workflow_config)
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_bad_name_format(self):
        """Ensure the function complains if the name doesn't match the expected Emu output format
        (so sample name can't be inferred)"""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "barcode02_RB02_emu.tsv"
        error_msg = "File name barcode02_RB02_emu.tsv does not conform to the expected format " \
                    "([SAMPLE]_[BARCODE]_rel-abundance.tsv). Sample name could not be extracted."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            summarize_emu.report_species_per_barcode(sample_path, self.workflow_config)

    def test_bad_sample_name(self):
        """Ensure the function complains if the name doesn't match the expected sample name format
        (so sample name can't be inferred)"""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "barcode2_RB02_rel-abundance.tsv"
        error_msg = "File name barcode2_RB02_rel-abundance.tsv does not conform to " \
                     "the expected format " \
                    "([SAMPLE]_[BARCODE]_rel-abundance.tsv). Sample name could not be extracted."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            summarize_emu.report_species_per_barcode(sample_path, self.workflow_config)

    def test_fail_wrong_abundance(self):
        """Ensure a relative abundance that does not sum to 1 (suggesting a corrupted file)
        is caught."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "barcode03_RB03_rel-abundance.tsv"
        error_msg = "Relative abundance does not sum to 1. " \
                    "This suggests the result file is broken (missing/extra lines)."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            summarize_emu.report_species_per_barcode(sample_path, self.workflow_config)


class TestMergeEmu(unittest.TestCase):
    def test_merge_identical_species(self):
        """Ensure dataframes with identical indexes can be merged."""
        barcode_1 = pd.DataFrame(data = {"abundance_from_all": [0.75, 0.2, 0.05],
                                         "estimated counts": [15.0, 4.0, 1.0],
                                         "medtages": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_1_header = ["barcode01_RB01"] * len(barcode_1.columns)
        barcode_1.columns = pd.MultiIndex.from_arrays([barcode_1_header,
                                                       barcode_1.columns])
        barcode_2 = pd.DataFrame(data = {"abundance_from_all": [0.8, 0.2, 0.00],
                                         "estimated counts": [16.0, 4.0, 0.0],
                                         "medtages": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_2_header = ["barcode02_RB02"] * len(barcode_2.columns)
        barcode_2.columns = pd.MultiIndex.from_arrays([barcode_2_header,
                                                       barcode_2.columns])
        expected_values = [[0.75, 15.0, "", 0.8, 16.0, ""],
                           [0.2, 4.0, "", 0.2, 4.0, ""],
                           [0.05, 1.0, "", 0.00, 0.0, ""]]
        expected_index = pd.Index(data = ["Placeholderia fakeorum",
                                          "Placeholderia bielefeldensis",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["barcode01_RB01",
                                                       "barcode01_RB01",
                                                       "barcode01_RB01",
                                                       "barcode02_RB02",
                                                       "barcode02_RB02",
                                                       "barcode02_RB02"],
                                                      ["abundance_from_all", "estimated counts",
                                                       "medtages",
                                                       "abundance_from_all", "estimated counts",
                                                       "medtages"]])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_emu([barcode_1, barcode_2])
        print(expected_merged)
        print(test_merged)
        pd.testing.assert_frame_equal(expected_merged, test_merged)

    def test_order_by_sample(self):
        """Ensure order of sample columns is consistent."""
        barcode_1 = pd.DataFrame(data = {"abundance_from_all": [0.75, 0.2, 0.05],
                                         "estimated counts": [15.0, 4.0, 1.0],
                                         "medtages": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_1_header = ["barcode01_RB01"] * len(barcode_1.columns)
        barcode_1.columns = pd.MultiIndex.from_arrays([barcode_1_header,
                                                       barcode_1.columns])
        barcode_2 = pd.DataFrame(data = {"abundance_from_all": [0.8, 0.2, 0.00],
                                         "estimated counts": [16.0, 4.0, 0.0],
                                         "medtages": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_2_header = ["barcode02_RB02"] * len(barcode_2.columns)
        barcode_2.columns = pd.MultiIndex.from_arrays([barcode_2_header,
                                                       barcode_2.columns])
        expected_values = [[0.75, 15.0, "", 0.8, 16.0, ""],
                           [0.2, 4.0, "", 0.2, 4.0, ""],
                           [0.05, 1.0, "", 0.00, 0.0, ""]]
        expected_index = pd.Index(data = ["Placeholderia fakeorum",
                                          "Placeholderia bielefeldensis",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["barcode01_RB01",
                                                       "barcode01_RB01",
                                                       "barcode01_RB01",
                                                       "barcode02_RB02",
                                                       "barcode02_RB02",
                                                       "barcode02_RB02"],
                                                      ["abundance_from_all", "estimated counts",
                                                       "medtages",
                                                       "abundance_from_all", "estimated counts",
                                                       "medtages"]])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_emu([barcode_2, barcode_1])
        print(expected_merged)
        print(test_merged)
        pd.testing.assert_frame_equal(expected_merged, test_merged)

    def test_merge_different_species(self):
        """Ensure dataframes with different indexes can be merged."""
        barcode_1 = pd.DataFrame(data = {"abundance_from_all": [0.75, 0.2, 0.05],
                                         "estimated counts": [15.0, 4.0, 1.0],
                                         "medtages": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_1_header = ["barcode01_RB01"] * len(barcode_1.columns)
        barcode_1.columns = pd.MultiIndex.from_arrays([barcode_1_header,
                                                       barcode_1.columns])
        barcode_2 = pd.DataFrame(data = {"abundance_from_all": [0.8, 0.2, 0.00],
                                         "estimated counts": [16.0, 4.0, 0.0],
                                         "medtages": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia testfacei",
                                                          "unassigned"], name = "species"))
        barcode_2_header = ["barcode02_RB02"] * len(barcode_2.columns)
        barcode_2.columns = pd.MultiIndex.from_arrays([barcode_2_header,
                                                       barcode_2.columns])
        expected_values = [[0.2, 4.0, "", np.nan, np.nan, np.nan],
                           [0.75, 15.0, "", 0.8, 16.0, ""],
                           [np.nan, np.nan, np.nan, 0.2, 4.0, ""],
                           [0.05, 1.0, "", 0.00, 0.0, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "Placeholderia testfacei",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["barcode01_RB01",
                                                       "barcode01_RB01",
                                                       "barcode01_RB01",
                                                       "barcode02_RB02",
                                                       "barcode02_RB02",
                                                       "barcode02_RB02"],
                                                      ["abundance_from_all", "estimated counts",
                                                       "medtages",
                                                       "abundance_from_all", "estimated counts",
                                                       "medtages"]])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_emu([barcode_1, barcode_2])
        print(expected_merged)
        print(test_merged)
        pd.testing.assert_frame_equal(expected_merged, test_merged)


    def test_multi_merge(self):
        """Merge more than two dataframes."""
        barcode_1 = pd.DataFrame(data = {"abundance_from_all": [0.75, 0.2, 0.05],
                                         "estimated counts": [15.0, 4.0, 1.0],
                                         "medtages": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_1_header = ["barcode01_RB01"] * len(barcode_1.columns)
        barcode_1.columns = pd.MultiIndex.from_arrays([barcode_1_header,
                                                       barcode_1.columns])
        barcode_2 = pd.DataFrame(data = {"abundance_from_all": [0.8, 0.2, 0.00],
                                         "estimated counts": [16.0, 4.0, 0.0],
                                         "medtages": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia testfacei",
                                                          "unassigned"], name = "species"))
        barcode_2_header = ["barcode02_RB02"] * len(barcode_2.columns)
        barcode_2.columns = pd.MultiIndex.from_arrays([barcode_2_header,
                                                       barcode_2.columns])
        barcode_3 = pd.DataFrame(data = {"abundance_from_all": [0.75, 0.2, 0.05],
                                         "estimated counts": [15.0, 4.0, 1.0],
                                         "medtages": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_3_header = ["barcode03_RB03"] * len(barcode_3.columns)
        barcode_3.columns = pd.MultiIndex.from_arrays([barcode_3_header,
                                                       barcode_3.columns])
        expected_values = [[0.2, 4.0, "", np.nan, np.nan, np.nan, 0.2, 4.0, ""],
                           [0.75, 15.0, "", 0.8, 16.0, "", 0.75, 15.0, ""],
                           [np.nan, np.nan, np.nan,  0.2, 4.0, "", np.nan, np.nan, np.nan,],
                           [0.05, 1.0, "", 0.00, 0.0, "", 0.05, 1.0, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "Placeholderia testfacei",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["barcode01_RB01",
                                                       "barcode01_RB01",
                                                       "barcode01_RB01",
                                                       "barcode02_RB02",
                                                       "barcode02_RB02",
                                                       "barcode02_RB02",
                                                       "barcode03_RB03",
                                                       "barcode03_RB03",
                                                       "barcode03_RB03"],
                                                      ["abundance_from_all", "estimated counts",
                                                       "medtages",
                                                       "abundance_from_all", "estimated counts",
                                                       "medtages",
                                                       "abundance_from_all", "estimated counts",
                                                       "medtages"
                                                       ]])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_emu([barcode_1, barcode_2, barcode_3])
        print(expected_merged)
        print(test_merged)
        pd.testing.assert_frame_equal(expected_merged, test_merged)


class TestMergeEmuDir(unittest.TestCase):
    workflow_config = {"sample_number_settings": {"sample_number_format": r"barcode\d{2}",
                                                  "positive_control": {},
                                                  "negative_control": ""},
                       "barcode_format": "RB[0-9]{2}"}

    def test_success_merge(self):
        """Test if multiple files are merged successfully."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "success_merge_dir"
        expected_values = [[0.2, 4.0, "", np.nan, np.nan, np.nan],
                           [0.75, 15.0, "", 0.8, 16.0, ""],
                           [np.nan, np.nan, np.nan, 0.2, 4.0, ""],
                           [0.05, 1.0, "", 0.00, 0.0, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "Placeholderia testfacei",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["barcode01_RB01",
                                                       "barcode01_RB01",
                                                       "barcode01_RB01",
                                                       "barcode02_RB02",
                                                       "barcode02_RB02",
                                                       "barcode02_RB02"],
                                                      ["abundance_from_all", "estimated counts",
                                                       "medtages",
                                                       "abundance_from_all", "estimated counts",
                                                       "medtages"]])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path,
                                                         active_config = self.workflow_config)
        print(expected_merged)
        print(test_merged)
        pd.testing.assert_frame_equal(expected_merged, test_merged)

    def test_handle_single_sample(self):
        """Test if a single file is parsed and returned properly."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "single_sample"
        expected_results = pd.DataFrame(data = {"abundance_from_all": [0.75, 0.2, 0.05],
                                                "estimated counts": [15.0, 4.0, 1.0],
                                                "medtages": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia fakeorum",
                                                                 "Placeholderia bielefeldensis",
                                                                 "unassigned"], name = "species"))
        barcode_header = ["barcode01_RB01"] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([barcode_header,
                                                              expected_results.columns])
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path,
                                                         active_config = self.workflow_config)
        pd.testing.assert_frame_equal(expected_results, test_merged)

    def test_handle_broken_file(self):
        """Ensure that an invalid file is handled properly (log error and return fake empty df)"""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "one_broken"
        expected_values = [[0.2, 4.0, "", np.nan, np.nan, np.nan],
                           [0.75, 15.0, "",np.nan, np.nan, np.nan],
                           [0.05, 1.0, "", np.nan, np.nan, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["barcode01_RB01",
                                                       "barcode01_RB01",
                                                       "barcode01_RB01",
                                                       "barcode03_RB03",
                                                       "barcode03_RB03",
                                                       "barcode03_RB03"],
                                                      ["abundance_from_all", "estimated counts",
                                                       "medtages",
                                                       "abundance_from_all", "estimated counts",
                                                       "medtages"]])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        with self.assertLogs("summarize_emu") as logged:
            log_msg = "ERROR:summarize_emu:Error in " \
                       f"{sample_path / 'barcode03_RB03_rel-abundance.tsv'}:\n" \
                      "Relative abundance does not sum to 1. " \
                      "This suggests the result file is broken (missing/extra lines).\n" \
                      "Empty results will be added to the merged summary."
            test_merged = summarize_emu.merge_all_in_emu_dir(sample_path,
                                                             active_config = self.workflow_config)
        print("\n".join(logged.output))
        assert log_msg in logged.output
        pd.testing.assert_frame_equal(expected_merged, test_merged)

    def test_merge_different_format_with_controls(self):
        """Ensure samples with different formats and controls are handled properly."""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "positive_control": {},
                                                      "negative_control": "NegK"},
                       "barcode_format": "RB[0-9]{2}"}
        sample_path = pathlib.Path(__file__).parent / "data"/"summarize_emu"/"merge_different_format"

        expected_values = [[0.75, 15.0, "", 0.75, 15.0, ""],
                           [0.2, 4.0, "", 0.2, 4.0, ""],
                           [0.05, 1.0, "", 0.05, 1.0, ""]]
        expected_index = pd.Index(data = ["Placeholderia fakeorum",
                                          "Placeholderia bielefeldensis",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["F99123456-0_RB01",
                                                       "F99123456-0_RB01",
                                                       "F99123456-0_RB01",
                                                       "NegK_RB02",
                                                       "NegK_RB02",
                                                       "NegK_RB02"],
                                                      ["abundance_from_all", "estimated counts",
                                                       "medtages",
                                                       "abundance_from_all", "estimated counts",
                                                       "medtages"]])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path,
                                                         active_config = workflow_config)
        print(expected_merged)
        pd.testing.assert_frame_equal(expected_merged, test_merged)

    def test_fail_no_files(self):
        """Ensure the function fails if no files matching the format are found."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "blank_dir"
        error_msg = f"No Emu reports found in {sample_path}."
        with pytest.raises(FileNotFoundError, match=re.escape(error_msg)):
            summarize_emu.merge_all_in_emu_dir(sample_path, active_config = self.workflow_config)


class TestWriteToSheets(unittest.TestCase):
    @classmethod
    def tearDownClass(cls) -> None:
        # remove test sheet
        (pathlib.Path(__file__).parent / "data" / "summarize_emu"
         / "test_results_sheet.xlsx").unlink()

    def test_write_success(self):
        expected_values = [[0.2, 4.0, "", np.nan, np.nan, np.nan],
                           [0.75, 15.0, "", 0.8, 16.0, ""],
                           [np.nan, np.nan, np.nan, 0.2, 4.0, ""],
                           [0.05, 1.0, "", 0.00, 0.0, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "Placeholderia testfacei",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["barcode01_RB01",
                                                       "barcode01_RB01",
                                                       "barcode01_RB01",
                                                       "barcode02_RB02",
                                                       "barcode02_RB02",
                                                       "barcode02_RB02"],
                                                      ["abundance_from_all", "estimated counts",
                                                       "medtages",
                                                       "abundance_from_all", "estimated counts",
                                                       "medtages"]])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        expected_abundance_values = [[0.2, np.nan],
                                     [0.75, 0.8],
                                     [np.nan, 0.2],
                                     [0.05, 0.00]]
        expected_abundance_cols = pd.MultiIndex.from_arrays([["barcode01_RB01",
                                                              "barcode02_RB02"],
                                                             ["abundance_from_all",
                                                              "abundance_from_all"]])
        expected_abundance = pd.DataFrame(data = expected_abundance_values, index = expected_index,
                                          columns = expected_abundance_cols)
        expected_count_values = [[4.0, np.nan],
                                 [15.0, 16.0],
                                 [np.nan, 4.0],
                                 [1.0, 0.00]]
        expected_count_cols = pd.MultiIndex.from_arrays([["barcode01_RB01",
                                                          "barcode02_RB02"],
                                                         ["estimated counts",
                                                          "estimated counts"]])
        expected_count = pd.DataFrame(data = expected_count_values, index = expected_index,
                                      columns = expected_count_cols)
        test_sheet = (pathlib.Path(__file__).parent / "data" / "summarize_emu"
         / "test_results_sheet.xlsx")
        summarize_emu.write_to_sheets(expected_merged, test_sheet)
        test_merged = pd.read_excel(test_sheet, sheet_name = "overview", index_col = 0,
                                    header = [0, 1])
        test_merged.loc[["Placeholderia bielefeldensis",
                         "Placeholderia fakeorum",
                         "unassigned"], pd.IndexSlice[["barcode01_RB01"],
                                                      ["medtages"]]] = ""
        test_merged.loc[["Placeholderia fakeorum",
                         "Placeholderia testfacei",
                         "unassigned"], pd.IndexSlice[["barcode02_RB02"],
                                                      ["medtages"]]] = ""
        test_abundance = pd.read_excel(test_sheet, sheet_name = "abundance", index_col = 0,
                                       header = [0, 1])
        test_count = pd.read_excel(test_sheet, sheet_name = "count", index_col = 0,
                                   header = [0, 1])
        pd.testing.assert_frame_equal(test_merged, expected_merged)
        pd.testing.assert_frame_equal(test_abundance, expected_abundance)
        pd.testing.assert_frame_equal(test_count, expected_count)

