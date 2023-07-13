import pathlib
import re
import unittest

import pandas as pd
import pytest

import summarize_emu


class TestExtractCounts(unittest.TestCase):
    def test_get_counts_success(self):
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "barcode01_rel-abundance.tsv"
        expected_results = pd.DataFrame(data={"species": ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"],
                                              "abundance_from_all": [0.75, 0.2, 0.05],
                                              "estimated counts": [15.0, 4.0, 1.0],
                                              "medtages": ["", "", ""]})
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

