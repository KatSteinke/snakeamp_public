import pathlib
import re
import unittest

import numpy as np
import pandas as pd
import pytest

import summarize_emu


class TestExtractCounts(unittest.TestCase):
    def test_get_counts_success(self):
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "barcode01_rel-abundance.tsv"
        expected_results = pd.DataFrame(data={"abundance_from_all": [0.75, 0.2, 0.05],
                                              "estimated counts": [15.0, 4.0, 1.0],
                                              "medtages": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia fakeorum",
                                                                 "Placeholderia bielefeldensis",
                                                                 "unassigned"], name = "species"))
        barcode_header = ["barcode01"] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([barcode_header,
                                                              expected_results.columns])
        test_results = summarize_emu.report_species_per_barcode(sample_path)
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_bad_name_format(self):
        """Ensure the function complains if the name doesn't match the expected Emu output format
        (so sample name can't be inferred)"""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "barcode02_emu.tsv"
        error_msg = "File name barcode02_emu.tsv does not conform to the expected format " \
                    "([SAMPLE]_rel-abundance.tsv). Sample name could not be extracted."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            summarize_emu.report_species_per_barcode(sample_path)

    def test_fail_wrong_abundance(self):
        """Ensure a relative abundance that does not sum to 1 (suggesting a corrupted file)
        is caught."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "barcode03_rel-abundance.tsv"
        error_msg = "Relative abundance does not sum to 1. " \
                    "This suggests the result file is broken (missing/extra lines)."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            summarize_emu.report_species_per_barcode(sample_path)


class TestMergeEmu(unittest.TestCase):
    def test_merge_identical_species(self):
        """Ensure dataframes with identical indexes can be merged."""
        barcode_1 = pd.DataFrame(data = {"abundance_from_all": [0.75, 0.2, 0.05],
                                         "estimated counts": [15.0, 4.0, 1.0],
                                         "medtages": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_1_header = ["barcode01"] * len(barcode_1.columns)
        barcode_1.columns = pd.MultiIndex.from_arrays([barcode_1_header,
                                                       barcode_1.columns])
        barcode_2 = pd.DataFrame(data = {"abundance_from_all": [0.8, 0.2, 0.00],
                                         "estimated counts": [16.0, 4.0, 0.0],
                                         "medtages": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_2_header = ["barcode02"] * len(barcode_2.columns)
        barcode_2.columns = pd.MultiIndex.from_arrays([barcode_2_header,
                                                       barcode_2.columns])
        expected_values = [[0.75, 15.0, "", 0.8, 16.0, ""],
                           [0.2, 4.0, "", 0.2, 4.0, ""],
                           [0.05, 1.0, "", 0.00, 0.0, ""]]
        expected_index = pd.Index(data = ["Placeholderia fakeorum",
                                          "Placeholderia bielefeldensis",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["barcode01", "barcode01", "barcode01",
                                                       "barcode02", "barcode02", "barcode02"],
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

    def test_merge_different_species(self):
        """Ensure dataframes with different indexes can be merged."""
        barcode_1 = pd.DataFrame(data = {"abundance_from_all": [0.75, 0.2, 0.05],
                                         "estimated counts": [15.0, 4.0, 1.0],
                                         "medtages": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_1_header = ["barcode01"] * len(barcode_1.columns)
        barcode_1.columns = pd.MultiIndex.from_arrays([barcode_1_header,
                                                       barcode_1.columns])
        barcode_2 = pd.DataFrame(data = {"abundance_from_all": [0.8, 0.2, 0.00],
                                         "estimated counts": [16.0, 4.0, 0.0],
                                         "medtages": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia testfacei",
                                                          "unassigned"], name = "species"))
        barcode_2_header = ["barcode02"] * len(barcode_2.columns)
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
        expected_columns = pd.MultiIndex.from_arrays([["barcode01", "barcode01", "barcode01",
                                                       "barcode02", "barcode02", "barcode02"],
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
        barcode_1_header = ["barcode01"] * len(barcode_1.columns)
        barcode_1.columns = pd.MultiIndex.from_arrays([barcode_1_header,
                                                       barcode_1.columns])
        barcode_2 = pd.DataFrame(data = {"abundance_from_all": [0.8, 0.2, 0.00],
                                         "estimated counts": [16.0, 4.0, 0.0],
                                         "medtages": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia testfacei",
                                                          "unassigned"], name = "species"))
        barcode_2_header = ["barcode02"] * len(barcode_2.columns)
        barcode_2.columns = pd.MultiIndex.from_arrays([barcode_2_header,
                                                       barcode_2.columns])
        barcode_3 = pd.DataFrame(data = {"abundance_from_all": [0.75, 0.2, 0.05],
                                         "estimated counts": [15.0, 4.0, 1.0],
                                         "medtages": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_3_header = ["barcode03"] * len(barcode_3.columns)
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
        expected_columns = pd.MultiIndex.from_arrays([["barcode01", "barcode01", "barcode01",
                                                       "barcode02", "barcode02", "barcode02",
                                                       "barcode03", "barcode03", "barcode03"],
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
