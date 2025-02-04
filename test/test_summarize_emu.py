import pathlib
import re
import unittest
from unittest import mock

import numpy as np
import pandas as pd
import pytest

import summarize_emu


class TestGetLISData(unittest.TestCase):
    workflow_config = {"sample_number_settings": {"sample_number_format":
                                                      r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                  "format_in_sheet":
                                                      r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                  "format_in_lis":
                                                      r'(?P<sample_type>[BDFT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                  "format_output": r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                  "positive_control": {"PosK": "Placeholderia"},
                                                  "negative_control": "NegK",
                                                  "sample_numbers_in": "number",
                                                  "sample_numbers_out": "letter",
                                                  "sample_numbers_output": "number",
                                                  "number_to_letter": {"70": "P", "30": "B",
                                                                       "10": "D", "50": "T"}
                                                  },
                       "barcode_format": "RB[0-9]{2}",
                       "lab_info_system": {"use_lis_features": True,
                                           "lis_report": (pathlib.Path(
                                               __file__).parent / "data" / "summarize_emu"
                                                          / "fake_mads_material.csv")}}
    lis_data = pd.read_csv(workflow_config["lab_info_system"]["lis_report"],
                           encoding = "latin1", dtype = {"modtaget": str,
                                                         "cprnr.": str})

    def test_fail_missing_number(self):
        """Fail if the sample number cannot be found in the LIS report."""
        sample_number = "1199123456-0"
        error_msg = ("Sample number F99123456 (original number: 1199123456-0)"
                     " not found in LIS report.")
        with pytest.raises(KeyError, match = re.escape(error_msg)):
            summarize_emu.get_lis_information(sample_number, self.lis_data, self.workflow_config)

    def test_get_data_success(self):
        """Correctly retrieve and convert data from LIS report."""
        expected_result = pd.DataFrame(data = {"patient": ["0000000000"],
                                               "prøvenr": ["F99123456"],
                                               "modtagedato": ["2021-01-02"],
                                               "prøvemateriale": ["Podning"],
                                               "anatomi": ["Svælg/tonsil"],
                                             "indikation": ["en eller anden lang tekst<Break/>"
                                                            "der ikke kan være på en linje i MADS"]})
        sample_number = "1199123456-0"
        test_result = summarize_emu.get_lis_information(sample_number, self.lis_data,
                                                        self.workflow_config)
        pd.testing.assert_frame_equal(expected_result, test_result)

    def test_get_blank_success(self):
        """Handle blank components in LIS report."""
        expected_result = pd.DataFrame(data = {"patient": ["0000000000"],
                                               "prøvenr": ["F99654321"],
                                               "modtagedato": ["2021-01-02"],
                                               "prøvemateriale": ["Spinalvæske"],
                                               "anatomi": [""],
                                             "indikation": ["en eller anden lang tekst<Break/>"
                                                            "der ikke kan være på en linje i MADS"]})
        sample_number = "1199654321-0"
        lis_data = pd.read_csv((pathlib.Path(__file__).parent / "data" / "summarize_emu"
                                / "fake_mads_material_blank.csv"),
                               encoding = "latin1", dtype = {"modtaget": str,
                                                             "cprnr.": str})
        test_result = summarize_emu.get_lis_information(sample_number, lis_data,
                                                        self.workflow_config)
        pd.testing.assert_frame_equal(expected_result, test_result)

    def test_handle_control(self):
        """Return blank results for controls."""
        expected_result = pd.DataFrame(data = {"patient": [""],
                                               "prøvenr": ["NegK"],
                                               "modtagedato": [""],
                                               "prøvemateriale": [""],
                                               "anatomi": [""],
                                             "indikation": [""]})
        sample_number = "NegK"
        test_result = summarize_emu.get_lis_information(sample_number, self.lis_data,
                                                        self.workflow_config)
        pd.testing.assert_frame_equal(expected_result, test_result)

        positive_expected = pd.DataFrame(data = {"patient": [""],
                                                 "prøvenr": ["PosK"],
                                                 "modtagedato": [""],
                                                 "prøvemateriale": [""],
                                                 "anatomi": [""],
                                                 "indikation": [""]})
        sample_number = "PosK"
        positive_test = summarize_emu.get_lis_information(sample_number, self.lis_data,
                                                          self.workflow_config)
        pd.testing.assert_frame_equal(positive_expected, positive_test)

    def test_handle_old_format(self):
        """Handle an older report without additional information on indication etc."""
        lis_data = pd.read_csv((pathlib.Path(__file__).parent / "data" / "summarize_emu"
                               / "test_mads_old.csv"),
                               encoding = "latin1", dtype = {"modtaget": str,
                                                             "cprnr.": str})
        expected_result = pd.DataFrame(data = {"patient": ["0000000000"],
                                               "prøvenr": ["F99123456"],
                                               "modtagedato": ["2021-01-02"],
                                               "prøvemateriale": ["Podning"],
                                               "anatomi": ["Svælg/tonsil"],
                                               "indikation": [""]})
        sample_number = "1199123456-0"
        test_result = summarize_emu.get_lis_information(sample_number, lis_data,
                                                        self.workflow_config)
        pd.testing.assert_frame_equal(expected_result, test_result)



class TestGetNameComponents(unittest.TestCase):
    workflow_config = {"sample_number_settings": {"sample_number_format": r"barcode\d{2}",
                                                  "positive_control": {},
                                                  "negative_control": ""},
                       "barcode_format": "RB[0-9]{2}",
                       "lab_info_system": {"use_lis_features": False}}

    def test_fail_bad_name(self):
        """Fail if the name doesn't conform to the expected format."""
        sample_name = "barcode02_RB02_emu.tsv"
        error_msg = "File name barcode02_RB02_emu.tsv does not conform to the expected format " \
                    "([RUN]_[SAMPLE]_[BARCODE]_rel-abundance.tsv)." \
                    " Sample name components could not be extracted."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            summarize_emu.extract_name_components(sample_name, self.workflow_config)

    def test_get_name(self):
        """Extract the components of the name where present."""
        sample_name = "kørsel0001-Y20231009-16S_barcode01_RB01_rel-abundance.tsv"
        expected_run = "kørsel0001-Y20231009-16S"
        expected_name = "barcode01"
        expected_barcode = "RB01"
        test_components = summarize_emu.extract_name_components(sample_name, self.workflow_config)
        assert expected_run == test_components.run_name
        assert expected_name == test_components.sample_name
        assert expected_barcode == test_components.barcode

    def test_get_short_name(self):
        """Extract the components of the name where present, using a shorter run number."""
        sample_name = "kørsel1-Y20231009-16S_barcode01_RB01_rel-abundance.tsv"
        expected_run = "kørsel1-Y20231009-16S"
        expected_name = "barcode01"
        expected_barcode = "RB01"
        test_components = summarize_emu.extract_name_components(sample_name, self.workflow_config)
        assert expected_run == test_components.run_name
        assert expected_name == test_components.sample_name
        assert expected_barcode == test_components.barcode

    def test_get_controls(self):
        """Get names for positive and negative controls."""
        test_config = {"input_names": pathlib.Path(__file__).parent / "data" / "input_names"
                                      / "input_da_old_lis.yaml",
                       "sample_number_settings": {"sample_number_format": r"barcode\d{2}",
                                                  "positive_control": {"PosK": ""},
                                                  "negative_control":
                                                      "NegK[a-zA-ZæøåÆØÅ0-9_-]*"},
                       "barcode_format": "RB[0-9]{2}",
                       "lab_info_system": {"use_lis_features": False}}
        posk_name_test = "kørsel0001-Y20231009-16S_PosK_RB01_rel-abundance.tsv"
        expected_posk_run = "kørsel0001-Y20231009-16S"
        expected_posk_name = "PosK"
        expected_posk_barcode = "RB01"
        test_components = summarize_emu.extract_name_components(posk_name_test,
                                                                test_config)
        assert expected_posk_run == test_components.run_name
        assert expected_posk_name == test_components.sample_name
        assert expected_posk_barcode == test_components.barcode
        negk_name_test = "kørsel0001-Y20231009-16S_NegK_RB01_rel-abundance.tsv"
        expected_negk_run = "kørsel0001-Y20231009-16S"
        expected_negk_name = "NegK"
        expected_negk_barcode = "RB01"
        test_components = summarize_emu.extract_name_components(negk_name_test,
                                                                test_config)
        assert expected_negk_run == test_components.run_name
        assert expected_negk_name == test_components.sample_name
        assert expected_negk_barcode == test_components.barcode
        
        # and now with one that overlaps with the sample format
        tricky_negk_name_test = "kørsel0001-Y20231009-16S_NegK__barcode01_RB01_rel-abundance.tsv"
        expected_tricky_negk_run = "kørsel0001-Y20231009-16S"
        expected_tricky_negk_name = "NegK__barcode01"
        expected_tricky_negk_barcode = "RB01"
        test_components = summarize_emu.extract_name_components(tricky_negk_name_test,
                                                                test_config)
        assert expected_tricky_negk_run == test_components.run_name
        assert expected_tricky_negk_name == test_components.sample_name
        assert expected_tricky_negk_barcode == test_components.barcode



class TestGetReadQC(unittest.TestCase):
    def test_success_all_qc(self):
        """Successfully extract Kraken and NanoStat results in a format that can be combined
        with Emu reports."""
        sample_dir = (pathlib.Path(__file__).parent / "data" / "summarize_emu"
                      / "result_base_dir" / "barcode01_RB01")
        expected_results = pd.DataFrame(data={"total_before_qc": [6176],
                                              "total_after_qc": [100],
                                              "human": [10]})
        test_results = summarize_emu.get_stats_from_sample_dir(sample_dir)
        pd.testing.assert_frame_equal(expected_results, test_results)


class TestExtractCounts(unittest.TestCase):
    workflow_config = {"sample_number_settings": {"sample_number_format": r"barcode\d{2}",
                                                  "positive_control": {},
                                                  "negative_control": ""},
                       "barcode_format": "RB[0-9]{2}",
                       "lab_info_system": {"use_lis_features": False}}
    base_dir = pathlib.Path(
        __file__).parent / "data" / "summarize_emu" / "result_base_dir"
    version_text = f"Version_{summarize_emu.__version__}"

    def test_bad_name_format(self):
        """Ensure the function complains if the name doesn't match the expected Emu output format
        (so sample name can't be inferred)"""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "barcode02_RB02_emu.tsv"
        error_msg = "File name barcode02_RB02_emu.tsv does not conform to the expected format " \
                    "([RUN]_[SAMPLE]_[BARCODE]_rel-abundance.tsv)." \
                    " Sample name components could not be extracted."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            summarize_emu.report_species_per_barcode(sample_path, self.base_dir,
                                                     self.workflow_config)

    def test_bad_sample_name(self):
        """Ensure the function complains if the name doesn't match the expected sample name format
        (so sample name can't be inferred)"""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "barcode2_RB02_rel-abundance.tsv"
        error_msg = "File name barcode2_RB02_rel-abundance.tsv does not conform to " \
                    "the expected format " \
                    "([RUN]_[SAMPLE]_[BARCODE]_rel-abundance.tsv)." \
                    " Sample name components could not be extracted."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            summarize_emu.report_species_per_barcode(sample_path, self.base_dir,
                                                     self.workflow_config)

    def test_fail_wrong_abundance(self):
        """Ensure a relative abundance that does not sum to 100% (suggesting a corrupted file)
        is caught."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "RUN0001_barcode03_RB03_rel-abundance.tsv"
        error_msg = "Relative abundance does not sum to 100%. " \
                    "This suggests the result file is broken (missing/extra lines)."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            summarize_emu.report_species_per_barcode(sample_path, self.base_dir,
                                                     self.workflow_config)

    def test_fail_wrong_abundance_single_taxid(self):
        """Ensure a relative abundance that does not sum to 100% (suggesting a corrupted file)
        is caught for a single taxid that isn't 'unassigned'."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "RUN0001_barcode04_RB04_rel-abundance.tsv"
        error_msg = "Relative abundance does not sum to 100%. " \
                    "This suggests the result file is broken (missing/extra lines)."
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            summarize_emu.report_species_per_barcode(sample_path, self.base_dir,
                                                     self.workflow_config)

    def test_get_counts_success(self):
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "RUN0001_barcode01_RB01_rel-abundance.tsv"
        expected_results = pd.DataFrame(data = {"abundance": [20.00, 75.00,
                                                              5.00],
                                                "counts" : [4, 15, 1],
                                                "med": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia bielefeldensis",
                                                                 "Placeholderia fakeorum",
                                                                 "unassigned"],
                                                         name = "species"))
        expected_results = expected_results.astype({"counts": "Int64"})
        run_header = ["RUN0001"] * len(expected_results.columns)
        version_header = [self.version_text] * len(expected_results.columns)
        name_header = ["barcode01"] * len(expected_results.columns)
        barcode_header = ["RB01"] * len(expected_results.columns)
        pre_qc_header = [6176] * len(expected_results.columns)
        post_qc_header = [100] * len(expected_results.columns)
        human_header = [10] * len(expected_results.columns)
        phhv_header = [""] * len(expected_results.columns)
        note_header = [""] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([run_header,
                                                              version_header,
                                                              barcode_header,
                                                              name_header,
                                                              pre_qc_header,
                                                              post_qc_header,
                                                              human_header,
                                                              phhv_header,
                                                              note_header,
                                                              expected_results.columns],
                                                             names = ["run",
                                                                      "pipeline_version",
                                                                      "barcode",
                                                                      "prøvenummer",
                                                                      "total_before_qc",
                                                                      "total_after_qc",
                                                                      "human",
                                                                      "PhHV",
                                                                      "notes",
                                                                      None])
        test_results = summarize_emu.report_species_per_barcode(sample_path, self.base_dir,
                                                                self.workflow_config)
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_get_counts_success_short_name(self):
        """Handle a run name without zero padding."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "RUN1_barcode01_RB01_rel-abundance.tsv"
        expected_results = pd.DataFrame(data = {"abundance": [20.00, 75.00,
                                                              5.00],
                                                "counts": [4, 15, 1],
                                                "med": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia bielefeldensis",
                                                                 "Placeholderia fakeorum",
                                                                 "unassigned"], name = "species"))
        expected_results = expected_results.astype({"counts": "Int64"})
        run_header = ["RUN1"] * len(expected_results.columns)
        version_header = [self.version_text] * len(expected_results.columns)
        name_header = ["barcode01"] * len(expected_results.columns)
        barcode_header = ["RB01"] * len(expected_results.columns)
        pre_qc_header = [6176] * len(expected_results.columns)
        post_qc_header = [100] * len(expected_results.columns)
        human_header = [10] * len(expected_results.columns)
        phhv_header = [""] * len(expected_results.columns)
        note_header = [""] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([run_header,
                                                              version_header,
                                                              barcode_header,
                                                              name_header,
                                                              pre_qc_header,
                                                              post_qc_header,
                                                              human_header,
                                                              phhv_header,
                                                              note_header,
                                                              expected_results.columns],
                                                             names = ["run",
                                                                      "pipeline_version",
                                                                      "barcode",
                                                                      "prøvenummer",
                                                                      "total_before_qc",
                                                                      "total_after_qc",
                                                                      "human",
                                                                      "PhHV",
                                                                      "notes",
                                                                      None])
        test_results = summarize_emu.report_species_per_barcode(sample_path, self.base_dir,
                                                                self.workflow_config)
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_handle_failed_sample_success(self):
        """Successfully parse a fallback file created when Emu fails."""
        sample_path = (pathlib.Path(__file__).parent / "data" / "summarize_emu" / "fallback_test"
                       / "RUN0001_barcode01_RB01_rel-abundance.tsv")
        expected_results = pd.DataFrame(data = {"abundance": [100.0],
                                                "counts": [100],
                                                "med": [""]},
                                        index = pd.Index(data = ["unassigned"],
                                                         name = "species"))
        expected_results = expected_results.astype({"counts": "Int64"})
        run_header = ["RUN0001"] * len(expected_results.columns)
        version_header = [self.version_text] * len(expected_results.columns)
        name_header = ["barcode01"] * len(expected_results.columns)
        barcode_header = ["RB01"] * len(expected_results.columns)
        pre_qc_header = [6176] * len(expected_results.columns)
        post_qc_header = [100] * len(expected_results.columns)
        human_header = [10] * len(expected_results.columns)
        phhv_header = [""] * len(expected_results.columns)
        note_header = [""] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([run_header,
                                                              version_header,
                                                              barcode_header,
                                                              name_header,
                                                              pre_qc_header,
                                                              post_qc_header,
                                                              human_header,
                                                              phhv_header,
                                                              note_header,
                                                              expected_results.columns],
                                                             names = ["run",
                                                                      "pipeline_version",
                                                                      "barcode",
                                                                      "prøvenummer",
                                                                      "total_before_qc",
                                                                      "total_after_qc",
                                                                      "human",
                                                                      "PhHV",
                                                                      "notes",
                                                                      None])
        log_msg = "INFO:summarize_emu:All reads for sample barcode01 are unassigned."
        with self.assertLogs("summarize_emu") as logged:
            test_results = summarize_emu.report_species_per_barcode(sample_path, self.base_dir,
                                                                    self.workflow_config)
            assert log_msg in logged.output
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_handle_blank_sample_success(self):
        """Successfully parse a blank fallback file created when Emu fails."""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "positive_control": {},
                                                      "negative_control": ""},
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": False}}
        sample_path = (pathlib.Path(__file__).parent / "data" / "summarize_emu" / "merge_one_empty"
                       / "RUN0001_1199123456-1_RB01_rel-abundance.tsv")
        expected_results = pd.DataFrame(data = {"abundance": [0.00],
                                                "counts": [0],
                                                "med": [""]},
                                        index = pd.Index(data = ["unassigned"],
                                                         name = "species"))
        expected_results = expected_results.astype({"counts": "Int64"})
        run_header = ["RUN0001"] * len(expected_results.columns)
        version_header = [self.version_text] * len(expected_results.columns)
        name_header = ["1199123456-1"] * len(expected_results.columns)
        barcode_header = ["RB01"] * len(expected_results.columns)
        pre_qc_header = [6176] * len(expected_results.columns)
        post_qc_header = [100] * len(expected_results.columns)
        human_header = [100] * len(expected_results.columns)
        phhv_header = [""] * len(expected_results.columns)
        note_header = [""] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([run_header,
                                                              version_header,
                                                              barcode_header,
                                                              name_header,
                                                              pre_qc_header,
                                                              post_qc_header,
                                                              human_header,
                                                              phhv_header,
                                                              note_header,
                                                              expected_results.columns],
                                                             names = ["run",
                                                                      "pipeline_version",
                                                                      "barcode",
                                                                      "prøvenummer",
                                                                      "total_before_qc",
                                                                      "total_after_qc",
                                                                      "human",
                                                                      "PhHV",
                                                                      "notes",
                                                                      None])
        log_msg = "INFO:summarize_emu:All reads for sample 1199123456-1 are unassigned."
        with self.assertLogs("summarize_emu") as logged:
            test_results = summarize_emu.report_species_per_barcode(sample_path, self.base_dir,
                                                                    workflow_config)
            assert log_msg in logged.output
        pd.testing.assert_frame_equal(expected_results, test_results, check_dtype = False)

    def test_handle_duplicate_orgs_success(self):
        sample_path = (pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "duplicate_orgs"
                       / "RUN0001_barcode01_RB01_rel-abundance.tsv")
        expected_results = pd.DataFrame(data = {"abundance": [20.00, 75.00,
                                                              5.00],
                                                "counts":[4, 15, 1],
                                                "med": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia bielefeldensis",
                                                                 "Placeholderia fakeorum",
                                                                 "unassigned"],
                                                         name = "species"))
        expected_results = expected_results.astype({"counts": "Int64"})
        run_header = ["RUN0001"] * len(expected_results.columns)
        version_header = [self.version_text] * len(expected_results.columns)
        name_header = ["barcode01"] * len(expected_results.columns)
        barcode_header = ["RB01"] * len(expected_results.columns)
        pre_qc_header = [6176] * len(expected_results.columns)
        post_qc_header = [100] * len(expected_results.columns)
        human_header = [10] * len(expected_results.columns)
        phhv_header = [""] * len(expected_results.columns)
        note_header = [""] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([run_header,
                                                              version_header,
                                                              barcode_header,
                                                              name_header,
                                                              pre_qc_header,
                                                              post_qc_header,
                                                              human_header,
                                                              phhv_header,
                                                              note_header,
                                                              expected_results.columns],
                                                             names = ["run",
                                                                      "pipeline_version",
                                                                      "barcode",
                                                                      "prøvenummer",
                                                                       "total_before_qc",
                                                                      "total_after_qc",
                                                                      "human",
                                                                      "PhHV",
                                                                      "notes",
                                                                      None])
        test_results = summarize_emu.report_species_per_barcode(sample_path, self.base_dir,
                                                                self.workflow_config)
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_handle_isolate_number(self):
        """Handle formats including special characters"""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      # TODO: replace with 0-9
                                                      "positive_control": {},
                                                      "negative_control": ""},
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": False}}
        sample_path = (pathlib.Path(
            __file__).parent / "data" / "summarize_emu"
                       / "RUN0001_F99123456-0_RB01_rel-abundance.tsv")
        expected_results = pd.DataFrame(data = {"abundance": [20.00, 75.00,
                                                              5.00],
                                                "counts": [4, 15, 1],
                                                "med": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia bielefeldensis",
                                                                 "Placeholderia fakeorum",
                                                                 "unassigned"],
                                                         name = "species"))
        expected_results = expected_results.astype({"counts": "Int64"})
        run_header = ["RUN0001"] * len(expected_results.columns)
        version_header = [self.version_text] * len(expected_results.columns)
        name_header = ["F99123456-0"] * len(expected_results.columns)
        barcode_header = ["RB01"] * len(expected_results.columns)
        pre_qc_header = [6176] * len(expected_results.columns)
        post_qc_header = [100] * len(expected_results.columns)
        human_header = [10] * len(expected_results.columns)
        phhv_header = [""] * len(expected_results.columns)
        note_header = [""] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([run_header,
                                                              version_header,
                                                              barcode_header,
                                                              name_header,
                                                              pre_qc_header,
                                                              post_qc_header,
                                                              human_header,
                                                              phhv_header,
                                                              note_header,
                                                              expected_results.columns],
                                                             names = ["run",
                                                                      "pipeline_version",
                                                                      "barcode",
                                                                      "prøvenummer",
                                                                      "total_before_qc",
                                                                      "total_after_qc",
                                                                      "human",
                                                                      "PhHV",
                                                                      "notes",
                                                                      None])
        test_results = summarize_emu.report_species_per_barcode(sample_path, self.base_dir,
                                                                workflow_config)
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_get_data_from_lis(self):
        """Optionally add metadata from LIS"""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "format_in_sheet":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "format_in_lis":
                                                          r'(?P<sample_type>[BDFT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                      "format_output":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "positive_control": {},
                                                      "negative_control": "",
                                                      "sample_numbers_in": "letter",
                                                      "sample_numbers_out": "letter",
                                                      "sample_numbers_output": "letter",
                                                      "number_to_letter": {"70": "P", "30": "B",
                                                                           "10": "D", "50": "T"}
                                                      },
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": True,
                                               "lis_report": (pathlib.Path(
                                                   __file__).parent / "data" / "summarize_emu"
                                                              / "fake_mads_material.csv")}}
        sample_path = (pathlib.Path(
            __file__).parent / "data" / "summarize_emu"
                       / "RUN0001_F99123456-0_RB01_rel-abundance.tsv")
        expected_results = pd.DataFrame(data = {"abundance": [20.00, 75.00,
                                                              5.00],
                                                "counts": [4, 15, 1],
                                                "med": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia bielefeldensis",
                                                                 "Placeholderia fakeorum",
                                                                 "unassigned"],
                                                         name = "species"))
        expected_results = expected_results.astype({"counts": "Int64"})
        run_header = ["RUN0001"] * len(expected_results.columns)
        version_header = [self.version_text] * len(expected_results.columns)
        name_header = ["F99123456"] * len(expected_results.columns)
        barcode_header = ["RB01"] * len(expected_results.columns)
        date_header = ["2021-01-02"] * len(expected_results.columns)
        patient_header = ["0000000000"] * len(expected_results.columns)
        material_header = ["Podning"] * len(expected_results.columns)
        anatomy_header = ["Svælg/tonsil"] * len(expected_results.columns)
        indication_header = ["en eller anden lang tekst<Break/>"
                             "der ikke kan være på en linje i MADS"] * len(expected_results.columns)
        pre_qc_header = [6176] * len(expected_results.columns)
        post_qc_header = [100] * len(expected_results.columns)
        human_header = [10] * len(expected_results.columns)
        phhv_header = [""] * len(expected_results.columns)
        note_header = [""] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([run_header,
                                                              version_header,
                                                              barcode_header,
                                                              name_header,
                                                              date_header,
                                                              patient_header,
                                                              material_header,
                                                              anatomy_header,
                                                              indication_header,
                                                              pre_qc_header,
                                                              post_qc_header,
                                                              human_header,
                                                              phhv_header,
                                                              note_header,
                                                              expected_results.columns],
                                                             names = ["run",
                                                                      "pipeline_version",
                                                                      "barcode",
                                                                      "prøvenummer",
                                                                      "modtagedato",
                                                                      "patient",
                                                                      "prøvemateriale",
                                                                      "anatomi",
                                                                      "indikation",
                                                                      "total_before_qc",
                                                                      "total_after_qc",
                                                                      "human",
                                                                      "PhHV",
                                                                      "notes",
                                                                      None]
                                                             )
        test_results = summarize_emu.report_species_per_barcode(sample_path, self.base_dir,
                                                                workflow_config)
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_translate_number(self):
        """Translate the sample number if needed"""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "format_in_sheet":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "format_in_lis":
                                                          r'(?P<sample_type>[BDFT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                      "format_output":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "positive_control": {},
                                                      "negative_control": "",
                                                      "sample_numbers_in": "number",
                                                      "sample_numbers_out": "letter",
                                                      "sample_numbers_output": "number",
                                                      "number_to_letter": {"70": "P", "30": "B",
                                                                           "10": "D", "50": "T"}
                                                      },
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": True,
                                               "lis_report": (pathlib.Path(
                                                   __file__).parent / "data" / "summarize_emu"
                                                              / "fake_mads_material.csv")}}
        sample_path = (pathlib.Path(
            __file__).parent / "data" / "summarize_emu"
                       / "RUN0001_1199123456-0_RB01_rel-abundance.tsv")

        expected_results = pd.DataFrame(data = {"abundance": [20.00, 75.00,
                                                              5.00],
                                                "counts":[4, 15, 1],
                                                "med": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia bielefeldensis",
                                                                 "Placeholderia fakeorum",
                                                                 "unassigned"],
                                                         name = "species"))
        expected_results = expected_results.astype({"counts": "Int64"})
        run_header = ["RUN0001"] * len(expected_results.columns)
        version_header = [self.version_text] * len(expected_results.columns)
        name_header = ["F99123456"] * len(expected_results.columns)
        barcode_header = ["RB01"] * len(expected_results.columns)
        date_header = ["2021-01-02"] * len(expected_results.columns)
        patient_header = ["0000000000"] * len(expected_results.columns)
        material_header = ["Podning"] * len(expected_results.columns)
        anatomy_header = ["Svælg/tonsil"] * len(expected_results.columns)
        indication_header = ["en eller anden lang tekst<Break/>"
                             "der ikke kan være på en linje i MADS"] * len(expected_results.columns)
        pre_qc_header = [6176] * len(expected_results.columns)
        post_qc_header = [100] * len(expected_results.columns)
        human_header = [10] * len(expected_results.columns)
        phhv_header = [""] * len(expected_results.columns)
        note_header = [""] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([run_header,
                                                              version_header,
                                                              barcode_header,
                                                              name_header,
                                                              date_header,
                                                              patient_header,
                                                              material_header,
                                                              anatomy_header,
                                                              indication_header,
                                                              pre_qc_header,
                                                              post_qc_header,
                                                              human_header,
                                                              phhv_header,
                                                              note_header,
                                                              expected_results.columns],
                                                             names = ["run",
                                                                      "pipeline_version",
                                                                      "barcode",
                                                                      "prøvenummer",
                                                                      "modtagedato",
                                                                      "patient",
                                                                      "prøvemateriale",
                                                                      "anatomi",
                                                                      "indikation",
                                                                       "total_before_qc",
                                                                      "total_after_qc",
                                                                      "human",
                                                                      "PhHV",
                                                                      "notes",
                                                                      None]
                                                             )
        test_results = summarize_emu.report_species_per_barcode(sample_path, self.base_dir,
                                                                workflow_config)
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_handle_control_material(self):
        """Insert blank "material" for controls"""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "format_in_sheet":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "format_in_lis":
                                                          r'(?P<sample_type>[BDFT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                      "format_output":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "positive_control": {},
                                                      "negative_control": "NegK",
                                                      "sample_numbers_in": "number",
                                                      "sample_numbers_out": "letter",
                                                      "sample_numbers_output": "number",
                                                      "number_to_letter": {"70": "P", "30": "B",
                                                                           "10": "D", "50": "T"}
                                                      },
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": True,
                                               "lis_report": (pathlib.Path(
                                                   __file__).parent / "data" / "summarize_emu"
                                                              / "fake_mads_material.csv")}}
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "RUN0001_NegK_RB02_rel-abundance.tsv"
        expected_results = pd.DataFrame(data = {"abundance": [20.00, 75.00,
                                                              5.00],
                                                "counts": [4, 15, 1],
                                                "med": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia bielefeldensis",
                                                                 "Placeholderia fakeorum",
                                                                 "unassigned"], name = "species"))
        expected_results = expected_results.astype({"counts": "Int64"})

        run_header = ["RUN0001"] * len(expected_results.columns)
        version_header = [self.version_text] * len(expected_results.columns)
        name_header = ["NegK"] * len(expected_results.columns)
        barcode_header = ["RB02"] * len(expected_results.columns)
        date_header = [""] * len(expected_results.columns)
        patient_header = [""] * len(expected_results.columns)
        material_header = [""] * len(expected_results.columns)
        anatomy_header = [""] * len(expected_results.columns)
        indication_header = [""] * len(expected_results.columns)
        pre_qc_header = [6176] * len(expected_results.columns)
        post_qc_header = [100] * len(expected_results.columns)
        human_header = [10] * len(expected_results.columns)
        phhv_header = [""] * len(expected_results.columns)
        note_header = [""] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([run_header,
                                                              version_header,
                                                              barcode_header,
                                                              name_header,
                                                              date_header,
                                                              patient_header,
                                                              material_header,
                                                              anatomy_header,
                                                              indication_header,
                                                              pre_qc_header,
                                                              post_qc_header,
                                                              human_header,
                                                              phhv_header,
                                                              note_header,
                                                              expected_results.columns],
                                                             names = ["run",
                                                                      "pipeline_version",
                                                                      "barcode",
                                                                      "prøvenummer",
                                                                      "modtagedato",
                                                                      "patient",
                                                                      "prøvemateriale",
                                                                      "anatomi",
                                                                      "indikation",
                                                                      "total_before_qc",
                                                                      "total_after_qc",
                                                                      "human",
                                                                      "PhHV",
                                                                      "notes",
                                                                      None])
        test_results = summarize_emu.report_species_per_barcode(sample_path, self.base_dir,
                                                                workflow_config)
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_handle_negative_control(self):
        """Handle negative controls"""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      # TODO: replace with 0-9
                                                      "positive_control": {},
                                                      "negative_control": "NegK"},
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": False}}
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "RUN0001_NegK_RB02_rel-abundance.tsv"
        expected_results = pd.DataFrame(data = {"abundance": [20.00, 75.00,
                                                              5.00],
                                                "counts": [4, 15, 1],
                                                "med": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia bielefeldensis",
                                                                 "Placeholderia fakeorum",
                                                                 "unassigned"], name = "species"))
        expected_results = expected_results.astype({"counts": "Int64"})

        run_header = ["RUN0001"] * len(expected_results.columns)
        version_header = [self.version_text] * len(expected_results.columns)
        name_header = ["NegK"] * len(expected_results.columns)
        barcode_header = ["RB02"] * len(expected_results.columns)
        pre_qc_header = [6176] * len(expected_results.columns)
        post_qc_header = [100] * len(expected_results.columns)
        human_header = [10] * len(expected_results.columns)
        phhv_header = [""] * len(expected_results.columns)
        note_header = [""] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([run_header,
                                                              version_header,
                                                              barcode_header,
                                                              name_header,
                                                              pre_qc_header,
                                                              post_qc_header,
                                                              human_header,
                                                              phhv_header,
                                                              note_header,
                                                              expected_results.columns],
                                                             names = ["run",
                                                                      "pipeline_version",
                                                                      "barcode",
                                                                      "prøvenummer",
                                                                      "total_before_qc",
                                                                      "total_after_qc",
                                                                      "human",
                                                                      "PhHV",
                                                                      "notes",
                                                                      None]
                                                             )
        test_results = summarize_emu.report_species_per_barcode(sample_path, self.base_dir,
                                                                workflow_config)
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_handle_positive_control(self):
        """Handle positive controls"""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      # TODO: replace with 0-9
                                                      "positive_control": {"PosK": "Placeholderia"},
                                                      "negative_control": ""},
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": False}}
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "RUN0001_PosK_RB03_rel-abundance.tsv"
        expected_results = pd.DataFrame(data = {"abundance": [20.00, 75.00,
                                                              5.00],
                                                "counts": [4, 15, 1],
                                                "med": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia bielefeldensis",
                                                                 "Placeholderia fakeorum",
                                                                 "unassigned"], name = "species"))
        expected_results = expected_results.astype({"counts": "Int64"})

        run_header = ["RUN0001"] * len(expected_results.columns)
        version_header = [self.version_text] * len(expected_results.columns)
        name_header = ["PosK"] * len(expected_results.columns)
        barcode_header = ["RB03"] * len(expected_results.columns)
        pre_qc_header = [6176] * len(expected_results.columns)
        post_qc_header = [100] * len(expected_results.columns)
        human_header = [10] * len(expected_results.columns)
        phhv_header = [""] * len(expected_results.columns)
        note_header = [""] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([run_header,
                                                              version_header,
                                                              barcode_header,
                                                              name_header,
                                                              pre_qc_header,
                                                              post_qc_header,
                                                              human_header,
                                                              phhv_header,
                                                              note_header,
                                                              expected_results.columns],
                                                             names = ["run",
                                                                      "pipeline_version",
                                                                      "barcode",
                                                                      "prøvenummer",
                                                                      "total_before_qc",
                                                                      "total_after_qc",
                                                                      "human",
                                                                      "PhHV",
                                                                      "notes",
                                                                      None])
        test_results = summarize_emu.report_species_per_barcode(sample_path, self.base_dir,
                                                                workflow_config)
        pd.testing.assert_frame_equal(expected_results, test_results)


class TestMergeEmu(unittest.TestCase):
    def test_merge_identical_species(self):
        """Ensure dataframes with identical indexes can be merged."""
        barcode_1 = pd.DataFrame(data = {"abundance": [75.00, 20.00, 5.00],
                                         "counts": [15, 4, 1],
                                         "med": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_1 = barcode_1.astype({"counts": "Int64"})

        name_1_header = ["barcode01"] * len(barcode_1.columns)
        barcode_1_header = ["RB01"] * len(barcode_1.columns)
        phhv_header_1 = [""] * len(barcode_1.columns)

        barcode_1.columns = pd.MultiIndex.from_arrays([barcode_1_header,
                                                       name_1_header,
                                                       phhv_header_1,
                                                       barcode_1.columns])
        barcode_2 = pd.DataFrame(data = {"abundance": [80.00, 20.00, 0.00],
                                         "counts": [16, 4, 0],
                                         "med": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_2 = barcode_2.astype({"counts": "Int64"})

        name_2_header = ["barcode02"] * len(barcode_2.columns)
        barcode_2_header = ["RB02"] * len(barcode_2.columns)
        phhv_header_2 = [""] * len(barcode_2.columns)

        barcode_2.columns = pd.MultiIndex.from_arrays([barcode_2_header,
                                                       name_2_header,
                                                       phhv_header_2,
                                                       barcode_2.columns])
        expected_values = [[20.00, 4, "", 20.00, 4, ""],
                           [75.00, 15, "", 80.00, 16, ""],
                           [5.00, 1, "", 0.00, 0, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02"],
                                                      ["barcode01",
                                                       "barcode01",
                                                       "barcode01",
                                                       "barcode02",
                                                       "barcode02",
                                                       "barcode02"],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_emu([barcode_1, barcode_2])
        print(expected_merged)
        print(test_merged)
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_order_by_sample(self):
        """Ensure order of sample columns is consistent."""
        barcode_1 = pd.DataFrame(data = {"abundance": [75.00, 20.00, 5.00],
                                         "counts": [15, 4, 1],
                                         "med": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_1 = barcode_1.astype({"counts": "Int64"})

        name_1_header = ["barcode01"] * len(barcode_1.columns)
        barcode_1_header = ["RB01"] * len(barcode_1.columns)
        phhv_header_1 = [""] * len(barcode_1.columns)
        barcode_1.columns = pd.MultiIndex.from_arrays([barcode_1_header, name_1_header,
                                                       phhv_header_1,
                                                       barcode_1.columns])
        barcode_2 = pd.DataFrame(data = {"abundance": [80.00, 20.00, 0.00],
                                         "counts": [16, 4, 0],
                                         "med": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_2 = barcode_2.astype({"counts": "Int64"})

        name_2_header = ["barcode02"] * len(barcode_2.columns)
        barcode_2_header = ["RB02"] * len(barcode_2.columns)
        phhv_header_2 = [""] * len(barcode_2.columns)
        barcode_2.columns = pd.MultiIndex.from_arrays([barcode_2_header, name_2_header,
                                                       phhv_header_2,
                                                       barcode_2.columns])
        expected_values = [[20.00, 4, "", 20.00, 4, ""],
                           [75.00, 15, "", 80.00, 16, ""],
                           [5.00, 1, "", 0.00, 0, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02"],
                                                      ["barcode01",
                                                       "barcode01",
                                                       "barcode01",
                                                       "barcode02",
                                                       "barcode02",
                                                       "barcode02"],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_emu([barcode_2, barcode_1])
        print(expected_merged)
        print(test_merged)
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_merge_different_species(self):
        """Ensure dataframes with different indexes can be merged."""
        barcode_1 = pd.DataFrame(data = {"abundance": [75.00, 20.00, 5.00],
                                         "counts": [15, 4, 1],
                                         "med": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_1 = barcode_1.astype({"counts": "Int64"})

        name_1_header = ["barcode01"] * len(barcode_1.columns)
        barcode_1_header = ["RB01"] * len(barcode_1.columns)
        phhv_header_1 = [""] * len(barcode_1.columns)
        barcode_1.columns = pd.MultiIndex.from_arrays([barcode_1_header, name_1_header,
                                                       phhv_header_1,
                                                       barcode_1.columns])
        barcode_2 = pd.DataFrame(data = {"abundance": [80.00, 20.00, 0.00],
                                         "counts": [16, 4, 0],
                                         "med": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia testfacei",
                                                          "unassigned"], name = "species"))
        barcode_2 = barcode_2.astype({"counts": "Int64"})

        name_2_header = ["barcode02"] * len(barcode_2.columns)
        barcode_2_header = ["RB02"] * len(barcode_2.columns)
        phhv_header_2 = [""] * len(barcode_2.columns)
        barcode_2.columns = pd.MultiIndex.from_arrays([barcode_2_header, name_2_header,
                                                       phhv_header_2,
                                                       barcode_2.columns])
        expected_values = [[20.00, 4, "", np.nan, np.nan, np.nan],
                           [75.00, 15, "", 80.00, 16, ""],
                           [np.nan, np.nan, np.nan, 20.00, 4, ""],
                           [5.00, 1, "", 0.00, 0, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "Placeholderia testfacei",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02"],
                                                      ["barcode01",
                                                       "barcode01",
                                                       "barcode01",
                                                       "barcode02",
                                                       "barcode02",
                                                       "barcode02"],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_emu([barcode_1, barcode_2])
        print(expected_merged)
        print(test_merged)
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_multi_merge(self):
        """Merge more than two dataframes."""
        barcode_1 = pd.DataFrame(data = {"abundance": [75.00, 20.00, 5.00],
                                         "counts": [15, 4, 1],
                                         "med": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_1 = barcode_1.astype({"counts": "Int64"})

        name_1_header = ["barcode01"] * len(barcode_1.columns)
        barcode_1_header = ["RB01"] * len(barcode_1.columns)
        phhv_header_1 = [""] * len(barcode_1.columns)
        barcode_1.columns = pd.MultiIndex.from_arrays([barcode_1_header, name_1_header,
                                                       phhv_header_1,
                                                       barcode_1.columns])
        barcode_2 = pd.DataFrame(data = {"abundance": [80.00, 20.00, 0.00],
                                         "counts": [16, 4, 0],
                                         "med": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia testfacei",
                                                          "unassigned"], name = "species"))
        barcode_2 = barcode_2.astype({"counts": "Int64"})

        name_2_header = ["barcode02"] * len(barcode_2.columns)
        barcode_2_header = ["RB02"] * len(barcode_2.columns)
        phhv_header_2 = [""] * len(barcode_2.columns)
        barcode_2.columns = pd.MultiIndex.from_arrays([barcode_2_header, name_2_header,
                                                       phhv_header_2,
                                                       barcode_2.columns])
        barcode_3 = pd.DataFrame(data = {"abundance": [75.00, 20.00, 5.00],
                                         "counts": [15, 4, 1],
                                         "med": ["", "", ""]},
                                 index = pd.Index(data = ["Placeholderia fakeorum",
                                                          "Placeholderia bielefeldensis",
                                                          "unassigned"], name = "species"))
        barcode_3 = barcode_3.astype({"counts": "Int64"})

        name_3_header = ["barcode03"] * len(barcode_3.columns)
        barcode_3_header = ["RB03"] * len(barcode_3.columns)
        phhv_header_3 = [""] * len(barcode_3.columns)
        barcode_3.columns = pd.MultiIndex.from_arrays([barcode_3_header, name_3_header,
                                                       phhv_header_3,
                                                       barcode_3.columns])
        expected_values = [[20.00, 4, "", np.nan, np.nan, np.nan, 20.00, 4, ""],
                           [75.00, 15, "", 80.00, 16, "", 75.00, 15, ""],
                           [np.nan, np.nan, np.nan, 20.00, 4, "", np.nan, np.nan, np.nan],
                           [5.00, 1, "", 0.00, 0, "", 5.00, 1, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "Placeholderia testfacei",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02",
                                                       "RB03",
                                                       "RB03",
                                                       "RB03"],
                                                      ["barcode01",
                                                       "barcode01",
                                                       "barcode01",
                                                       "barcode02",
                                                       "barcode02",
                                                       "barcode02",
                                                       "barcode03",
                                                       "barcode03",
                                                       "barcode03"],
                                                      ["", "", "",
                                                       "", "", "",
                                                       "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"
                                                       ]])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_emu([barcode_1, barcode_2, barcode_3])
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)


class TestSortColumns(unittest.TestCase):
    species_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                     "Placeholderia fakeorum",
                                     "Placeholderia testfacei",
                                     "unassigned"], name = "species")

    def test_sort_no_controls(self):
        """Sort a report with no controls."""
        workflow_config = {"sample_number_settings": {"sample_number_format": r"sample\d{2}",
                                                      "positive_control": {},
                                                      "negative_control": ""},
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": False}}
        expected_values = [[20.00, 4, "", np.nan, np.nan, np.nan],
                           [75.00, 15, "", 80.00, 16, ""],
                           [np.nan, np.nan, np.nan, 20.00, 4, ""],
                           [5.00, 1, "", 0.00, 0, ""]]
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02"],
                                                      ["sample02",
                                                       "sample02",
                                                       "sample02",
                                                       "sample01",
                                                       "sample01",
                                                       "sample01"],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "PhHV",
                                                              None])
        expected_sorted = pd.DataFrame(data = expected_values, index = self.species_index,
                                       columns = expected_columns)
        unsorted_values = [[np.nan, np.nan, np.nan, 20.00, 4, ""],
                           [80.00, 16, "", 75.00, 15, ""],
                           [20.00, 4, "", np.nan, np.nan, np.nan],
                           [0.00, 0, "", 5.00, 1, ""]]
        unsorted_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      ["RB02",
                                                       "RB02",
                                                       "RB02",
                                                       "RB01",
                                                       "RB01",
                                                       "RB01"],
                                                      ["sample01",
                                                       "sample01",
                                                       "sample01",
                                                       "sample02",
                                                       "sample02",
                                                       "sample02"],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "PhHV",
                                                              None])
        test_df = pd.DataFrame(data = unsorted_values, index = self.species_index,
                               columns = unsorted_columns)
        test_sort = summarize_emu.sort_report_samples(test_df, workflow_config)
        pd.testing.assert_frame_equal(test_sort, expected_sorted)

    def test_sort_one_control_type(self):
        """Sort a report with only one type of controls."""
        workflow_config = {"sample_number_settings":
                               {"sample_number_format":
                                    r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                "positive_control": {},
                                "negative_control": "NegK"},
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": False}}
        species_counts = [[20.00, 4, "", 20.00, 4, "", 20.00, 4, ""],
                          [75.00, 15, "", 75.00, 15, "", 75.00, 15, ""],
                          [5.00, 1, "", 5.00, 1, "", 5.00, 1, ""],
                          [0.00, 0, "", 0.00, 0, "", 0.00, 0, ""]]
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      ["RB02",
                                                       "RB02",
                                                       "RB02",
                                                       "RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB03",
                                                       "RB03",
                                                       "RB03"],
                                                      ["NegK",
                                                       "NegK",
                                                       "NegK",
                                                       "F99123456-1",
                                                       "F99123456-1",
                                                       "F99123456-1",
                                                       "F99123456-0",
                                                       "F99123456-0",
                                                       "F99123456-0"],
                                                      ["", "", "", "", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run", "barcode",
                                                              "prøvenummer", "PhHV", None])
        expected_sorted = pd.DataFrame(data = species_counts, index = self.species_index,
                                       columns = expected_columns)
        unsorted_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02",
                                                       "RB03",
                                                       "RB03",
                                                       "RB03"],
                                                      ["F99123456-1",
                                                       "F99123456-1",
                                                       "F99123456-1",
                                                       "NegK",
                                                       "NegK",
                                                       "NegK",
                                                       "F99123456-0",
                                                       "F99123456-0",
                                                       "F99123456-0"],
                                                      ["", "", "", "", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run", "barcode",
                                                              "prøvenummer", "PhHV", None])
        test_df = pd.DataFrame(data = species_counts, index = self.species_index,
                               columns = unsorted_columns)
        test_sort = summarize_emu.sort_report_samples(test_df, workflow_config)
        pd.testing.assert_frame_equal(test_sort, expected_sorted)

    def test_sort_with_duplicates(self):
        """Sort a report with duplicated sample numbers."""
        workflow_config = {"sample_number_settings":
                               {"sample_number_format":
                                    r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                "positive_control": {},
                                "negative_control": "NegK"},
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": False}}
        species_counts_sorted = [[20.00, 4, "", 20.00, 4, "", 25.00, 5, "", 20.00, 4, ""],
                                 [75.00, 15, "", 75.00, 15, "", 75.00, 15, "", 75.00, 15, ""],
                                 [5.00, 1, "", 5.00, 1, "", 0.00, 0, "", 5.00, 1, ""],
                                 [0.00, 0, "", 0.00, 0, "", 0.00, 0, "", 0.00, 0, ""]]
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      ["RB02",
                                                       "RB02",
                                                       "RB02",
                                                       "RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB03",
                                                       "RB03",
                                                       "RB03",
                                                       "RB04",
                                                       "RB04",
                                                       "RB04"],
                                                      ["NegK",
                                                       "NegK",
                                                       "NegK",
                                                       "F99123456-1",
                                                       "F99123456-1",
                                                       "F99123456-1",
                                                       "F99123456-0",
                                                       "F99123456-0",
                                                       "F99123456-0",
                                                       "F99123456-1",
                                                       "F99123456-1",
                                                       "F99123456-1"],
                                                      ["", "", "",
                                                       "", "", "",
                                                       "", "", "",
                                                       "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run", "barcode",
                                                              "prøvenummer", "PhHV", None])
        expected_sorted = pd.DataFrame(data = species_counts_sorted, index = self.species_index,
                                       columns = expected_columns)
        unsorted_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02",
                                                       "RB03",
                                                       "RB03",
                                                       "RB03",
                                                       "RB04",
                                                       "RB04",
                                                       "RB04"],
                                                      ["F99123456-1",
                                                       "F99123456-1",
                                                       "F99123456-1",
                                                       "NegK",
                                                       "NegK",
                                                       "NegK",
                                                       "F99123456-0",
                                                       "F99123456-0",
                                                       "F99123456-0",
                                                       "F99123456-1",
                                                       "F99123456-1",
                                                       "F99123456-1"],
                                                      ["", "", "",
                                                       "", "", "",
                                                       "", "", "",
                                                       "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run", "barcode",
                                                              "prøvenummer", "PhHV", None])
        test_df = pd.DataFrame(data = species_counts_sorted, index = self.species_index,
                               columns = unsorted_columns)
        test_sort = summarize_emu.sort_report_samples(test_df, workflow_config)
        pd.testing.assert_frame_equal(test_sort, expected_sorted)

    def test_sort_with_controls(self):
        """Sort a report with both positive and negative controls."""
        workflow_config = {"sample_number_settings":
                               {"sample_number_format":
                                    r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                "positive_control": {"PosK": "Placeholderia"},
                                "negative_control": "NegK"},
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": False}}
        species_counts = [[20.00, 4, "", 20.00, 4, "", 20.00, 4, ""],
                          [75.00, 15, "", 75.00, 15, "", 75.00, 15, ""],
                          [5.00, 1, "", 5.00, 1, "", 5.00, 1, ""],
                          [0.00, 0, "", 0.00, 0, "", 0.00, 0, ""]]
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      ["RB02",
                                                       "RB02",
                                                       "RB02",
                                                       "RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB03",
                                                       "RB03",
                                                       "RB03"],
                                                      ["NegK",
                                                       "NegK",
                                                       "NegK",
                                                       "PosK",
                                                       "PosK",
                                                       "PosK",
                                                       "F99123456-0",
                                                       "F99123456-0",
                                                       "F99123456-0"],
                                                      ["", "", "", "", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run", "barcode",
                                                              "prøvenummer", "PhHV", None])
        expected_sorted = pd.DataFrame(data = species_counts, index = self.species_index,
                                       columns = expected_columns)
        unsorted_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02",
                                                       "RB03",
                                                       "RB03",
                                                       "RB03"],
                                                      ["PosK",
                                                       "PosK",
                                                       "PosK",
                                                       "NegK",
                                                       "NegK",
                                                       "NegK",
                                                       "F99123456-0",
                                                       "F99123456-0",
                                                       "F99123456-0"],
                                                      ["", "", "", "", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run", "barcode",
                                                              "prøvenummer", "PhHV", None])
        test_df = pd.DataFrame(data = species_counts, index = self.species_index,
                               columns = unsorted_columns)
        test_sort = summarize_emu.sort_report_samples(test_df, workflow_config)
        pd.testing.assert_frame_equal(test_sort, expected_sorted)


class TestMergeEmuDir(unittest.TestCase):
    workflow_config = {"sample_number_settings": {"sample_number_format": r"barcode\d{2}",
                                                  "positive_control": {},
                                                  "negative_control": ""},
                       "barcode_format": "RB[0-9]{2}",
                       "lab_info_system": {"use_lis_features": False}}
    base_dir = pathlib.Path(__file__).parent / "data" / "summarize_emu" / "result_base_dir"
    version_text = f"Version_{summarize_emu.__version__}"

    def test_fail_missing_dir(self):
        """Complain informatively when the directory to merge doesn't exist."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "no_such_dir"
        error_msg = f"Emu report directory {sample_path} does not exist"
        with pytest.raises(FileNotFoundError, match = error_msg):
            summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                               active_config = self.workflow_config)

    def test_fail_not_a_dir_dir(self):
        """Complain informatively when pointed at a single file."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "barcode2_RB02_rel-abundance.tsv"
        error_msg = (f"{sample_path} is a single file. "
                     "Please specify the directory containing all Emu reports.")
        with pytest.raises(NotADirectoryError, match = error_msg):
            summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                               active_config = self.workflow_config)

    def test_success_merge(self):
        """Test if multiple files are merged successfully."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "success_merge_dir"
        expected_values = [[20.00, 4, "", np.nan, np.nan, np.nan],
                           [75.00, 15, "", 80.00, 16, ""],
                           [np.nan, np.nan, np.nan, 20.00, 4, ""],
                           [5.00, 1, "", 0.00, 0, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "Placeholderia testfacei",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02"],
                                                      ["barcode01",
                                                       "barcode01",
                                                       "barcode01",
                                                       "barcode02",
                                                       "barcode02",
                                                       "barcode02"],
                                                      [6176, 6176, 6176, 6176, 6176, 6176],
                                                      [100, 100, 100, 100, 100, 100],
                                                      [10, 10, 10, 10, 10, 10],
                                                      ["", "", "", "", "", ""],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                              "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                         active_config = self.workflow_config)
        print(expected_merged)
        print(test_merged)
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_handle_single_sample(self):
        """Test if a single file is parsed and returned properly."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "single_sample"
        expected_results = pd.DataFrame(data = {"abundance": [20.00, 75.00,
                                                              5.00],
                                                "counts": [4, 15,
                                                           1],
                                                "med": ["", "", ""]},
                                        index = pd.Index(data = ["Placeholderia bielefeldensis",
                                                                 "Placeholderia fakeorum",
                                                                 "unassigned"], name = "species"))
        expected_results = expected_results.astype({"counts": "Int64"})
        run_header = ["RUN0001"] * len(expected_results.columns)
        version_header = [self.version_text] * len(expected_results.columns)
        name_header = ["barcode01"] * len(expected_results.columns)
        barcode_header = ["RB01"] * len(expected_results.columns)
        pre_qc_header = [6176] * len(expected_results.columns)
        post_qc_header = [100] * len(expected_results.columns)
        human_header = [10] * len(expected_results.columns)
        phhv_header = [""] * len(expected_results.columns)
        note_header = [""] * len(expected_results.columns)
        expected_results.columns = pd.MultiIndex.from_arrays([run_header,
                                                              version_header,
                                                              barcode_header, name_header,
                                                              pre_qc_header,
                                                              post_qc_header,
                                                              human_header,
                                                              phhv_header,
                                                              note_header,
                                                              expected_results.columns],
                                                             names = ["run",
                                                                      "pipeline_version",
                                                                      "barcode",
                                                                      "prøvenummer",
                                                                      "total_before_qc",
                                                                      "total_after_qc",
                                                                      "human",
                                                                      "PhHV",
                                                                      "notes",
                                                                      None])
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                         active_config = self.workflow_config)
        pd.testing.assert_frame_equal(expected_results, test_merged)

    def test_handle_bad_read_data(self):
        """Handle a file for which both read QC data and emu report are broken."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "one_broken"
        base_dir = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "bad_qc"
        expected_values = [[20.00, 4, "", np.nan, np.nan, np.nan],
                           [75.00, 15, "", np.nan, np.nan, np.nan],
                           [5.00, 1, "", np.nan, np.nan, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB03",
                                                       "RB03",
                                                       "RB03"],
                                                      ["barcode01",
                                                       "barcode01",
                                                       "barcode01",
                                                       "barcode03",
                                                       "barcode03",
                                                       "barcode03"],
                                                      [6176, 6176, 6176, np.nan, np.nan, np.nan],
                                                      [100, 100, 100, np.nan, np.nan, np.nan],
                                                      [10, 10, 10, np.nan, np.nan, np.nan],
                                                      ["", "", "", "", "", ""],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                               "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        with self.assertLogs("summarize_emu") as logged:
            log_msg = "ERROR:summarize_emu:Error in " \
                      f"{sample_path / 'RUN0001_barcode03_RB03_rel-abundance.tsv'}:\n" \
                      "Relative abundance does not sum to 100%. " \
                      "This suggests the result file is broken (missing/extra lines).\n" \
                      "Empty results will be added to the merged summary."
            bad_output_msg = ("ERROR:summarize_emu:Error extracting QC data for barcode03:\n"
                              "'Kraken report is malformed. "
                              "Check that you are supplying a --report file.'\n"
                              "No QC data will be added.")
            test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, base_dir,
                                                             active_config = self.workflow_config)
        print("\n".join(logged.output))
        assert log_msg in logged.output
        assert bad_output_msg in logged.output
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False,
                                      check_column_type = False)

    def test_handle_missing_read_data(self):
        """Handle a file for which the emu report is broken and read data is missing."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "one_broken"
        base_dir =  pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "missing_qc"
        expected_values = [[20.00, 4, "", np.nan, np.nan, np.nan],
                           [75.00, 15, "", np.nan, np.nan, np.nan],
                           [5.00, 1, "", np.nan, np.nan, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB03",
                                                       "RB03",
                                                       "RB03"],
                                                      ["barcode01",
                                                       "barcode01",
                                                       "barcode01",
                                                       "barcode03",
                                                       "barcode03",
                                                       "barcode03"],
                                                      [6176, 6176, 6176, np.nan, np.nan, np.nan],
                                                      [100, 100, 100, np.nan, np.nan, np.nan],
                                                      [10, 10, 10, np.nan, np.nan, np.nan],
                                                      ["", "", "", "", "", ""],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                              "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        with self.assertLogs("summarize_emu") as logged:
            log_msg = "ERROR:summarize_emu:Error in " \
                      f"{sample_path / 'RUN0001_barcode03_RB03_rel-abundance.tsv'}:\n" \
                      "Relative abundance does not sum to 100%. " \
                      "This suggests the result file is broken (missing/extra lines).\n" \
                      "Empty results will be added to the merged summary."
            missing_data_msg = ("ERROR:summarize_emu:Error extracting QC data for barcode03:\n"
                                "Kraken report"
                                f" {base_dir / 'barcode03_RB03' / 'reads' / 'barcode03_RB03.kraken.tsv'}"
                                " not found.\n"
                                "No QC data will be added.")
            test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, base_dir,
                                                             active_config = self.workflow_config)
        assert log_msg in logged.output
        assert missing_data_msg in logged.output
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False,
                                      check_column_type = False)

    def test_handle_broken_file(self):
        """Ensure that an invalid file is handled properly (log error and return fake empty df)"""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "one_broken"
        expected_values = [[20.00, 4, "", np.nan, np.nan, np.nan],
                           [75.00, 15, "", np.nan, np.nan, np.nan],
                           [5.00, 1, "", np.nan, np.nan, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB03",
                                                       "RB03",
                                                       "RB03"],
                                                      ["barcode01",
                                                       "barcode01",
                                                       "barcode01",
                                                       "barcode03",
                                                       "barcode03",
                                                       "barcode03"],
                                                      [6176, 6176, 6176, 6176, 6176, 6176],
                                                      [100, 100, 100, 100, 100, 100],
                                                      [10, 10, 10, 10, 10, 10],
                                                      ["", "", "", "", "", ""],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                              "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        with self.assertLogs("summarize_emu") as logged:
            log_msg = "ERROR:summarize_emu:Error in " \
                      f"{sample_path / 'RUN0001_barcode03_RB03_rel-abundance.tsv'}:\n" \
                      "Relative abundance does not sum to 100%. " \
                      "This suggests the result file is broken (missing/extra lines).\n" \
                      "Empty results will be added to the merged summary."
            test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                             active_config = self.workflow_config)
        assert log_msg in logged.output
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_warn_different_runs(self):
        """Warn when a directory contains data from multiple runs."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "multi_run"
        expected_values = [[20.00, 4, "", np.nan, np.nan, np.nan],
                           [75.00, 15, "", 80.00, 16, ""],
                           [np.nan, np.nan, np.nan, 20.00, 4, ""],
                           [5.00, 1, "", 0.00, 0, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "Placeholderia testfacei",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001", "RUN0001", "RUN0001",
                                                       "RUN0002", "RUN0002", "RUN0002"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02"],
                                                      ["barcode01",
                                                       "barcode01",
                                                       "barcode01",
                                                       "barcode02",
                                                       "barcode02",
                                                       "barcode02"],
                                                      [6176, 6176, 6176, 6176, 6176, 6176],
                                                      [100, 100, 100, 100, 100, 100],
                                                      [10, 10, 10, 10, 10, 10],
                                                      ["", "", "", "", "", ""],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                              "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        with self.assertLogs("summarize_emu") as logged:
            log_msg = ("WARNING:summarize_emu:Data appear to be from multiple runs "
                       "(['RUN0001', 'RUN0002']).")
            test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                             active_config = self.workflow_config)
        print("\n".join(logged.output))
        assert log_msg in logged.output
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_handle_broken_with_lis(self):
        """Ensure that an invalid file is handled properly when using LIS data"""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "format_in_sheet":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "format_in_lis":
                                                          r'(?P<sample_type>[BDFT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                      "format_output":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "positive_control": {},
                                                      "negative_control": "NegK",
                                                      "sample_numbers_in": "letter",
                                                      "sample_numbers_out": "letter",
                                                      "sample_numbers_output": "letter",
                                                      "number_to_letter": {"70": "P", "30": "B",
                                                                           "10": "D", "50": "T"}
                                                      },
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": True,
                                               "lis_report": (pathlib.Path(
                                                   __file__).parent / "data" / "summarize_emu"
                                                              / "fake_mads_material.csv")}}
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "one_broken_lis"
        expected_values = [[20.00, 4, "", np.nan, np.nan, np.nan],
                           [75.00, 15, "", np.nan, np.nan, np.nan],
                           [5.00, 1, "", np.nan, np.nan, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02"],
                                                      ["F99123456",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99654321",
                                                       "F99654321",
                                                       "F99654321"],
                                                      ["2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "",
                                                       "",
                                                       ""],
                                                      ["RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "",
                                                       "",
                                                       ""],
                                                      ["Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "",
                                                       "",
                                                       ""],
                                                      ["Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "",
                                                       "",
                                                       ""],
                                                      ["en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "",
                                                       "",
                                                       ""],
                                                      [6176, 6176, 6176, 6176, 6176, 6176],
                                                      [100, 100, 100, 100, 100, 100],
                                                      [10, 10, 10, 10, 10, 10],
                                                      ["", "", "", "", "", ""],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "modtagedato",
                                                              "patient",
                                                              "prøvemateriale",
                                                              "anatomi",
                                                              "indikation",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                              "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        with self.assertLogs("summarize_emu") as logged:
            log_msg = "ERROR:summarize_emu:Error in " \
                      f"{sample_path / 'RUN0001_F99654321-0_RB02_rel-abundance.tsv'}:\n" \
                      "Relative abundance does not sum to 100%. " \
                      "This suggests the result file is broken (missing/extra lines).\n" \
                      "Empty results will be added to the merged summary."
            test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                             active_config = workflow_config)
        print("\n".join(logged.output))
        assert log_msg in logged.output
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_merge_different_format_with_controls(self):
        """Ensure samples with different formats and controls are handled properly."""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "positive_control": {},
                                                      "negative_control": "NegK"},
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": False}}
        sample_path = (pathlib.Path(__file__).parent / "data" / "summarize_emu"
                       / "merge_different_format")

        expected_values = [[20.00, 4, "", 20.00, 4, ""],
                           [75.00, 15, "", 75.00, 15, ""],
                           [5.00, 1, "", 5.00, 1, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB02",
                                                       "RB02",
                                                       "RB02",
                                                       "RB01",
                                                       "RB01",
                                                       "RB01"],
                                                      ["NegK",
                                                       "NegK",
                                                       "NegK",
                                                       "F99123456-0",
                                                       "F99123456-0",
                                                       "F99123456-0"],
                                                      [6176, 6176, 6176, 6176, 6176, 6176],
                                                      [100, 100, 100, 100, 100, 100],
                                                      [10, 10, 10, 10, 10, 10],
                                                      ["", "", "", "", "", ""],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                               "PhHV",
                                                              "notes", None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                         active_config = workflow_config)
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_multi_sort_with_controls(self):
        """Ensure non-control samples are sorted by barcode."""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "positive_control": {},
                                                      "negative_control": "NegK"},
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": False}}
        sample_path = (pathlib.Path( __file__).parent / "data" / "summarize_emu"
                       / "merge_multi_sample")

        expected_values = [[20.00, 4, "", 20.00, 4, "", 20.00, 4, ""],
                           [75.00, 15, "", 75.00, 15, "", 75.00, 15, ""],
                           [5.00, 1, "", 5.00, 1, "", 5.00, 1, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB02",
                                                       "RB02",
                                                       "RB02",
                                                       "RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB03",
                                                       "RB03",
                                                       "RB03"],
                                                      ["NegK",
                                                       "NegK",
                                                       "NegK",
                                                       "F99123456-1",
                                                       "F99123456-1",
                                                       "F99123456-1",
                                                       "F99123456-0",
                                                       "F99123456-0",
                                                       "F99123456-0"],
                                                      [6176, 6176, 6176,
                                                       6176, 6176, 6176,
                                                       6176, 6176, 6176],
                                                      [100, 100, 100,
                                                       100, 100, 100,
                                                       100, 100, 100],
                                                      [10, 10, 10,
                                                       10, 10, 10,
                                                       10, 10, 10],
                                                      ["", "", "", "", "", "", "", "", ""],
                                                      ["", "", "", "", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                              "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                         active_config = workflow_config)
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_handle_blank(self):
        """Handle a blank sample file."""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "format_in_sheet":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "format_in_lis":
                                                          r'(?P<sample_type>[BDFT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                      "format_output":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "positive_control": {},
                                                      "negative_control": "NegK",
                                                      "sample_numbers_in": "number",
                                                      "sample_numbers_out": "letter",
                                                      "sample_numbers_output": "number",
                                                      "number_to_letter": {"70": "P", "30": "B",
                                                                           "10": "D", "50": "T"}
                                                      },
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": True,
                                               "lis_report": (pathlib.Path(
                                                   __file__).parent / "data" / "summarize_emu"
                                                              / "fake_mads_material.csv")}}
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "merge_one_empty"

        expected_values = [[20.00, 4, "", np.nan, np.nan, np.nan, 20.00, 4, ""],
                           [75.00, 15, "", np.nan, np.nan, np.nan,  75.00, 15, ""],
                           [5.00, 1, "", 0.00, 0.00, "", 5.00, 1, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB02",
                                                       "RB02",
                                                       "RB02",
                                                       "RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB03",
                                                       "RB03",
                                                       "RB03"],
                                                      ["NegK",
                                                       "NegK",
                                                       "NegK",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99123456"],
                                                      ["",
                                                       "",
                                                       "",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02"],
                                                      ["",
                                                       "",
                                                       "",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0"],
                                                      ["",
                                                       "",
                                                       "",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning"
                                                       ],
                                                      ["",
                                                       "",
                                                       "",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       ],
                                                      ["", "", "",
                                                          "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS"],
                                                      [6176, 6176, 6176,
                                                       6176, 6176, 6176,
                                                       6176, 6176, 6176],
                                                      [100, 100, 100,
                                                       100, 100, 100,
                                                       100, 100, 100],
                                                      [10, 10, 10,
                                                       100, 100, 100,
                                                       10, 10, 10],
                                                      ["", "", "", "", "", "", "", "", ""],
                                                      ["", "", "", "", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "modtagedato",
                                                              "patient",
                                                              "prøvemateriale",
                                                              "anatomi",
                                                              "indikation",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                              "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                         active_config = workflow_config)
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)


    def test_handle_multi_blank(self):
        """Handle multiple blank sample files."""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "format_in_sheet":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "format_in_lis":
                                                          r'(?P<sample_type>[BDFT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                      "format_output":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "positive_control": {},
                                                      "negative_control": "NegK",
                                                      "sample_numbers_in": "number",
                                                      "sample_numbers_out": "letter",
                                                      "sample_numbers_output": "number",
                                                      "number_to_letter": {"70": "P", "30": "B",
                                                                           "10": "D", "50": "T"}
                                                      },
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": True,
                                               "lis_report": (pathlib.Path(
                                                   __file__).parent / "data" / "summarize_emu"
                                                              / "fake_mads_material.csv")}}
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "merge_two_empty"

        expected_values = [[20.00, 4, "", np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                           [75.00, 15, "", np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                           [5.00, 1, "", 0.00, 0.00, "", 0.00, 0, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB02",
                                                       "RB02",
                                                       "RB02",
                                                       "RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB03",
                                                       "RB03",
                                                       "RB03"],
                                                      ["NegK",
                                                       "NegK",
                                                       "NegK",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99123456"],
                                                      ["",
                                                       "",
                                                       "",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02"],
                                                      ["",
                                                       "",
                                                       "",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0"],
                                                      ["",
                                                       "",
                                                       "",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning"
                                                       ],
                                                      ["",
                                                       "",
                                                       "",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       ],
                                                      ["", "", "",
                                                          "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS"],
                                                      [6176, 6176, 6176,
                                                       6176, 6176, 6176,
                                                       6176, 6176, 6176],
                                                      [100, 100, 100,
                                                       100, 100, 100,
                                                       100, 100, 100],
                                                      [10, 10, 10,
                                                       100, 100, 100,
                                                       10, 10, 10],
                                                      ["", "", "", "", "", "", "", "", ""],
                                                      ["", "", "", "", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "modtagedato",
                                                              "patient",
                                                              "prøvemateriale",
                                                              "anatomi",
                                                              "indikation",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                              "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                         active_config = workflow_config)
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_merge_and_get_material(self):
        """Get sample material for all samples."""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "format_in_sheet":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "format_in_lis":
                                                          r'(?P<sample_type>[BDFT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                      "format_output":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "positive_control": {},
                                                      "negative_control": "NegK",
                                                      "sample_numbers_in": "letter",
                                                      "sample_numbers_out": "letter",
                                                      "sample_numbers_output": "letter",
                                                      "number_to_letter": {"70": "P", "30": "B",
                                                                           "10": "D", "50": "T"}
                                                      },
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": True,
                                               "lis_report": (pathlib.Path(
                                                   __file__).parent / "data" / "summarize_emu"
                                                              / "fake_mads_material.csv")}}
        sample_path = (pathlib.Path(__file__).parent / "data" / "summarize_emu"
                       / "merge_different_format")

        expected_values = [[20.00, 4, "", 20.00, 4, ""],
                           [75.00, 15, "", 75.00, 15, ""],
                           [5.00, 1, "", 5.00, 1, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB02",
                                                       "RB02",
                                                       "RB02",
                                                       "RB01",
                                                       "RB01",
                                                       "RB01"],
                                                      ["NegK",
                                                       "NegK",
                                                       "NegK",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99123456"],
                                                      ["",
                                                       "",
                                                       "",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       ],
                                                      ["",
                                                       "",
                                                       "",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0"],
                                                      ["",
                                                       "",
                                                       "",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning"
                                                       ],
                                                      ["",
                                                       "",
                                                       "",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       ],
                                                      ["",
                                                       "",
                                                       "",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS"],
                                                      [6176, 6176, 6176, 6176, 6176, 6176],
                                                      [100, 100, 100, 100, 100, 100],
                                                      [10, 10, 10, 10, 10, 10],
                                                      ["", "", "", "", "", ""],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "modtagedato",
                                                              "patient",
                                                              "prøvemateriale",
                                                              "anatomi",
                                                              "indikation",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                              "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                         active_config = workflow_config)
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_get_different_patients(self):
        """Ensure samples from different patients are reported correctly."""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "format_in_sheet":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "format_in_lis":
                                                          r'(?P<sample_type>[BDFT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                      "format_output":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "positive_control": {},
                                                      "negative_control": "NegK",
                                                      "sample_numbers_in": "letter",
                                                      "sample_numbers_out": "letter",
                                                      "sample_numbers_output": "letter",
                                                      "number_to_letter": {"70": "P", "30": "B",
                                                                           "10": "D", "50": "T"}
                                                      },
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": True,
                                               "lis_report": (pathlib.Path(
                                                   __file__).parent / "data" / "summarize_emu"
                                                              / "fake_mads_patients.csv")}}
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "merge_blank_material"

        expected_values = [[20.00, 4, "", 20.00, 4, ""],
                           [75.00, 15, "", 75.00, 15, ""],
                           [5.00, 1, "", 5.00, 1, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02"],
                                                      ["F99123456",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99654321",
                                                       "F99654321",
                                                       "F99654321"],
                                                      ["2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       ],
                                                      ["RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_1",
                                                       "RUN0001_pt_1",
                                                       "RUN0001_pt_1"],
                                                      ["Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning"
                                                       ],
                                                      ["Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil"],
                                                      ["en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS"],
                                                      [6176, 6176, 6176, 6176, 6176, 6176],
                                                      [100, 100, 100, 100, 100, 100],
                                                      [10, 10, 10, 10, 10, 10],
                                                      ["", "", "", "", "", ""],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "modtagedato",
                                                              "patient",
                                                              "prøvemateriale",
                                                              "anatomi",
                                                              "indikation",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                              "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                         active_config = workflow_config)
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_get_different_patients_different_runs(self):
        """Ensure samples from different patients from different runs are reported correctly."""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "format_in_sheet":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "format_in_lis":
                                                          r'(?P<sample_type>[BDFT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                      "format_output":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "positive_control": {},
                                                      "negative_control": "NegK",
                                                      "sample_numbers_in": "letter",
                                                      "sample_numbers_out": "letter",
                                                      "sample_numbers_output": "letter",
                                                      "number_to_letter": {"70": "P", "30": "B",
                                                                           "10": "D", "50": "T"}
                                                      },
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": True,
                                               "lis_report": (pathlib.Path(
                                                   __file__).parent / "data" / "summarize_emu"
                                                              / "fake_mads_patients.csv")}}
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "multi_run_patients"

        expected_values = [[20.00, 4, "", 20.00, 4, ""],
                           [75.00, 15, "", 75.00, 15, ""],
                           [5.00, 1, "", 5.00, 1, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0002",
                                                       "RUN0002",
                                                       "RUN0002"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02"],
                                                      ["F99123456",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99654321",
                                                       "F99654321",
                                                       "F99654321"],
                                                      ["2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       ],
                                                      ["RUN0001_RUN0002_pt_0",
                                                       "RUN0001_RUN0002_pt_0",
                                                       "RUN0001_RUN0002_pt_0",
                                                       "RUN0001_RUN0002_pt_1",
                                                       "RUN0001_RUN0002_pt_1",
                                                       "RUN0001_RUN0002_pt_1"],
                                                      ["Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning"
                                                       ],
                                                      ["Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil"],
                                                      ["en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS"],
                                                      [6176, 6176, 6176, 6176, 6176, 6176],
                                                      [100, 100, 100, 100, 100, 100],
                                                      [10, 10, 10, 10, 10, 10],
                                                      ["", "", "", "", "", ""],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "modtagedato",
                                                              "patient",
                                                              "prøvemateriale",
                                                              "anatomi",
                                                              "indikation",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                              "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                         active_config = workflow_config)
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_get_same_patient(self):
        """Ensure samples from the same patient are reported correctly."""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "format_in_sheet":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "format_in_lis":
                                                          r'(?P<sample_type>[BDFT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                      "format_output":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "positive_control": {},
                                                      "negative_control": "NegK",
                                                      "sample_numbers_in": "letter",
                                                      "sample_numbers_out": "letter",
                                                      "sample_numbers_output": "letter",
                                                      "number_to_letter": {"70": "P", "30": "B",
                                                                           "10": "D", "50": "T"}
                                                      },
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": True,
                                               "lis_report": (pathlib.Path(
                                                   __file__).parent / "data" / "summarize_emu"
                                                              / "fake_mads_material.csv")}}
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "merge_blank_material"

        expected_values = [[20.00, 4, "", 20.00, 4, ""],
                           [75.00, 15, "", 75.00, 15, ""],
                           [5.00, 1, "", 5.00, 1, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001",
                                                       "RUN0001"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02"],
                                                      ["F99123456",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99654321",
                                                       "F99654321",
                                                       "F99654321"],
                                                      ["2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       ],
                                                      ["RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0",
                                                       "RUN0001_pt_0"],
                                                      ["Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning"
                                                       ],
                                                      ["Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil"],
                                                      ["en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS"],
                                                      [6176, 6176, 6176, 6176, 6176, 6176],
                                                      [100, 100, 100, 100, 100, 100],
                                                      [10, 10, 10, 10, 10, 10],
                                                      ["", "", "", "", "", ""],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "modtagedato",
                                                              "patient",
                                                              "prøvemateriale",
                                                              "anatomi",
                                                              "indikation",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                              "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                         active_config = workflow_config)
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_shorten_run_name(self):
        """Ensure run name is shortened properly."""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "format_in_sheet":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "format_in_lis":
                                                          r'(?P<sample_type>[BDFT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                      "format_output":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "positive_control": {},
                                                      "negative_control": "NegK",
                                                      "sample_numbers_in": "letter",
                                                      "sample_numbers_out": "letter",
                                                      "sample_numbers_output": "letter",
                                                      "number_to_letter": {"70": "P", "30": "B",
                                                                           "10": "D", "50": "T"}
                                                      },
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": True,
                                               "lis_report": (pathlib.Path(
                                                   __file__).parent / "data" / "summarize_emu"
                                                              / "fake_mads_material.csv")}}
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "merge_shorten_name"

        expected_values = [[20.00, 4, "", 20.00, 4, ""],
                           [75.00, 15, "", 75.00, 15, ""],
                           [5.00, 1, "", 5.00, 1, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["Y20990101_Run0001",
                                                       "Y20990101_Run0001",
                                                       "Y20990101_Run0001",
                                                       "Y20990101_Run0001",
                                                       "Y20990101_Run0001",
                                                       "Y20990101_Run0001"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02"],
                                                      ["F99123456",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99654321",
                                                       "F99654321",
                                                       "F99654321"],
                                                      ["2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       ],
                                                      ["Run0001_pt_0",
                                                       "Run0001_pt_0",
                                                       "Run0001_pt_0",
                                                       "Run0001_pt_0",
                                                       "Run0001_pt_0",
                                                       "Run0001_pt_0"],
                                                      ["Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning"
                                                       ],
                                                      ["Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil"],
                                                      ["en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS"],
                                                      [6176, 6176, 6176, 6176, 6176, 6176],
                                                      [100, 100, 100, 100, 100, 100],
                                                      [10, 10, 10, 10, 10, 10],
                                                      ["", "", "", "", "", ""],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "modtagedato",
                                                              "patient",
                                                              "prøvemateriale",
                                                              "anatomi",
                                                              "indikation",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                              "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                         active_config = workflow_config)
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_handle_new_format(self):
        """Ensure run names without zero padding are handled properly."""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "format_in_sheet":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "format_in_lis":
                                                          r'(?P<sample_type>[BDFT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                      "format_output":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "positive_control": {},
                                                      "negative_control": "NegK",
                                                      "sample_numbers_in": "letter",
                                                      "sample_numbers_out": "letter",
                                                      "sample_numbers_output": "letter",
                                                      "number_to_letter": {"70": "P", "30": "B",
                                                                           "10": "D", "50": "T"}
                                                      },
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": True,
                                               "lis_report": (pathlib.Path(
                                                   __file__).parent / "data" / "summarize_emu"
                                                              / "fake_mads_material.csv")}}
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "merge_short_run_name"

        expected_values = [[20.00, 4, "", 20.00, 4, ""],
                           [75.00, 15, "", 75.00, 15, ""],
                           [5.00, 1, "", 5.00, 1, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["Y20990101_Run001",
                                                       "Y20990101_Run001",
                                                       "Y20990101_Run001",
                                                       "Y20990101_Run001",
                                                       "Y20990101_Run001",
                                                       "Y20990101_Run001"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02"],
                                                      ["F99123456",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99654321",
                                                       "F99654321",
                                                       "F99654321"],
                                                      ["2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       ],
                                                      ["Run001_pt_0",
                                                       "Run001_pt_0",
                                                       "Run001_pt_0",
                                                       "Run001_pt_0",
                                                       "Run001_pt_0",
                                                       "Run001_pt_0"],
                                                      ["Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning"
                                                       ],
                                                      ["Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil"],
                                                      ["en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS"],
                                                      [6176, 6176, 6176, 6176, 6176, 6176],
                                                      [100, 100, 100, 100, 100, 100],
                                                      [10, 10, 10, 10, 10, 10],
                                                      ["", "", "", "", "", ""],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "modtagedato",
                                                              "patient",
                                                              "prøvemateriale",
                                                              "anatomi",
                                                              "indikation",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                              "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                         active_config = workflow_config)
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_handle_broken_run_name(self):
        """Ensure a run name not matching the expected format is handled."""
        workflow_config = {"sample_number_settings": {"sample_number_format":
                                                          r'([BDFT]|[135]0|11)([0-9]{8}|[0-9]{6})-\d',
                                                      "format_in_sheet":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "format_in_lis":
                                                          r'(?P<sample_type>[BDFT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                      "format_output":
                                                          r'(?P<sample_type>[BDFT]|[135]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)',
                                                      "positive_control": {},
                                                      "negative_control": "NegK",
                                                      "sample_numbers_in": "letter",
                                                      "sample_numbers_out": "letter",
                                                      "sample_numbers_output": "letter",
                                                      "number_to_letter": {"70": "P", "30": "B",
                                                                           "10": "D", "50": "T"}
                                                      },
                           "barcode_format": "RB[0-9]{2}",
                           "lab_info_system": {"use_lis_features": True,
                                               "lis_report": (pathlib.Path(
                                                   __file__).parent / "data" / "summarize_emu"
                                                              / "fake_mads_material.csv")}}
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "merge_broken_run_name"

        expected_values = [[20.00, 4, "", 20.00, 4, ""],
                           [75.00, 15, "", 75.00, 15, ""],
                           [5.00, 1, "", 5.00, 1, ""]]
        expected_index = pd.Index(data = ["Placeholderia bielefeldensis",
                                          "Placeholderia fakeorum",
                                          "unassigned"], name = "species")
        expected_columns = pd.MultiIndex.from_arrays([["Y20990101_0001",
                                                       "Y20990101_0001",
                                                       "Y20990101_0001",
                                                       "Y20990101_0001",
                                                       "Y20990101_0001",
                                                       "Y20990101_0001"],
                                                      [self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text,
                                                       self.version_text],
                                                      ["RB01",
                                                       "RB01",
                                                       "RB01",
                                                       "RB02",
                                                       "RB02",
                                                       "RB02"],
                                                      ["F99123456",
                                                       "F99123456",
                                                       "F99123456",
                                                       "F99654321",
                                                       "F99654321",
                                                       "F99654321"],
                                                      ["2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       "2021-01-02",
                                                       ],
                                                      ["RUNxxxx_pt_0",
                                                       "RUNxxxx_pt_0",
                                                       "RUNxxxx_pt_0",
                                                       "RUNxxxx_pt_0",
                                                       "RUNxxxx_pt_0",
                                                       "RUNxxxx_pt_0"],
                                                      ["Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning",
                                                       "Podning"
                                                       ],
                                                      ["Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil",
                                                       "Svælg/tonsil"],
                                                      ["en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS",
                                                       "en eller anden lang tekst<Break/>"
                                                       "der ikke kan være på en linje i MADS"],
                                                      [6176, 6176, 6176, 6176, 6176, 6176],
                                                      [100, 100, 100, 100, 100, 100],
                                                      [10, 10, 10, 10, 10, 10],
                                                      ["", "", "", "", "", ""],
                                                      ["", "", "", "", "", ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["run",
                                                              "pipeline_version",
                                                              "barcode",
                                                              "prøvenummer",
                                                              "modtagedato",
                                                              "patient",
                                                              "prøvemateriale",
                                                              "anatomi",
                                                              "indikation",
                                                              "total_before_qc",
                                                              "total_after_qc",
                                                              "human",
                                                              "PhHV",
                                                              "notes",
                                                              None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        test_merged = summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                                         active_config = workflow_config)
        pd.testing.assert_frame_equal(expected_merged, test_merged, check_dtype = False)

    def test_fail_no_files(self):
        """Ensure the function fails if no files matching the format are found."""
        sample_path = pathlib.Path(
            __file__).parent / "data" / "summarize_emu" / "blank_dir"
        error_msg = f"No Emu reports found in {sample_path}."
        with pytest.raises(FileNotFoundError, match = re.escape(error_msg)):
            summarize_emu.merge_all_in_emu_dir(sample_path, self.base_dir,
                                               active_config = self.workflow_config)


class TestWriteToSheets(unittest.TestCase):
    test_sheet = (pathlib.Path(__file__).parent / "data" / "summarize_emu"
                  / "test_results_sheet.xlsx")
    expected_values = [[20, 4, "", np.nan, np.nan, np.nan],
                       [75, 15, "", 80, 16, ""],
                       [np.nan, np.nan, np.nan, 20, 4, ""],
                       [5, 1, "", 0, 0, ""]]
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
                                                  ["abundance", "counts",
                                                   "med",
                                                   "abundance", "counts",
                                                   "med"]],
                                                 names = ["prøvenummer", None])
    expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                   columns = expected_columns)
    print(expected_merged.columns.names)

    @classmethod
    def tearDownClass(cls) -> None:
        # remove test sheet
        (pathlib.Path(__file__).parent / "data" / "summarize_emu"
         / "test_results_sheet.xlsx").unlink()
        (pathlib.Path(__file__).parent / "data" / "summarize_emu"
         / "test_results_mads.xlsx").unlink()
        (pathlib.Path(__file__).parent / "data" / "summarize_emu"
         / "test_results_mads_blank.xlsx").unlink()

    def test_check_all_sheets_present(self):
        """Ensure that all sheets are present."""
        expected_sheet_names = ["overview", "abundance", "count", "notes"]
        summarize_emu.write_to_sheets(self.expected_merged, self.test_sheet)
        with pd.ExcelFile(self.test_sheet) as test_sheet_data:
            assert test_sheet_data.sheet_names == expected_sheet_names

    def test_add_notes_sheet(self):
        """Check that the notes sheet (with sample numbers as index) is created correctly."""
        summarize_emu.write_to_sheets(self.expected_merged, self.test_sheet)
        expected_notes = pd.DataFrame(index = pd.Index(data = ["barcode01_RB01", "barcode02_RB02"],
                                                       name = "Prøvenummer"),
                                      columns = ["notes"])
        test_notes = pd.read_excel(self.test_sheet, sheet_name = "notes", index_col = 0)
        pd.testing.assert_frame_equal(expected_notes, test_notes, check_dtype = False)

    def test_write_success(self):
        expected_abundance_values = [[20.00, np.nan],
                                     [75.00, 80.00],
                                     [np.nan, 20.00],
                                     [5.00, 0.00]]
        expected_abundance_cols = pd.MultiIndex.from_arrays([["barcode01_RB01",
                                                              "barcode02_RB02"],
                                                             ["abundance",
                                                              "abundance"]],
                                                            names = ["prøvenummer", None])
        expected_abundance = pd.DataFrame(data = expected_abundance_values,
                                          index = self.expected_index,
                                          columns = expected_abundance_cols)
        expected_count_values = [[4, np.nan],
                                 [15, 16],
                                 [np.nan, 4],
                                 [1, 0]]
        expected_count_cols = pd.MultiIndex.from_arrays([["barcode01_RB01",
                                                          "barcode02_RB02"],
                                                         ["counts",
                                                          "counts"]],
                                                        names = ["prøvenummer", None])
        expected_count = pd.DataFrame(data = expected_count_values,
                                      index = self.expected_index,
                                      columns = expected_count_cols)

        summarize_emu.write_to_sheets(self.expected_merged, self.test_sheet)
        test_merged = pd.read_excel(self.test_sheet, sheet_name = "overview", index_col = 0,
                                    header = [0, 1])
        test_merged.loc[["Placeholderia bielefeldensis",
                         "Placeholderia fakeorum",
                         "unassigned"], pd.IndexSlice[["barcode01_RB01"],
        ["med"]]] = ""
        test_merged.loc[["Placeholderia fakeorum",
                         "Placeholderia testfacei",
                         "unassigned"], pd.IndexSlice[["barcode02_RB02"],
        ["med"]]] = ""
        test_abundance = pd.read_excel(self.test_sheet, sheet_name = "abundance", index_col = 0,
                                       header = [0, 1])
        test_count = pd.read_excel(self.test_sheet, sheet_name = "count", index_col = 0,
                                   header = [0, 1])
        pd.testing.assert_frame_equal(test_merged, self.expected_merged, check_dtype = False)
        pd.testing.assert_frame_equal(test_abundance, expected_abundance, check_dtype = False)
        pd.testing.assert_frame_equal(test_count, expected_count, check_dtype = False)

    def test_handle_extra_lines(self):
        """Handle extra lines in the multiindex."""
        expected_values = [[20.00, 4, "", np.nan, np.nan, np.nan],
                           [75.00, 15, "", 80.00, 16, ""],
                           [np.nan, np.nan, np.nan, 20.00, 4, ""],
                           [5.00, 1, "", 0.00, 0, ""]]
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
                                                      ["podning",
                                                       "podning",
                                                       "podning",
                                                       "Væv",
                                                       "Væv",
                                                       "Væv"],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["prøvenummer",
                                                              "prøvemateriale", None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        expected_abundance_values = [[20.00, np.nan],
                                     [75.00, 80.00],
                                     [np.nan, 20.00],
                                     [5.00, 0.00]]
        expected_abundance_cols = pd.MultiIndex.from_arrays([["barcode01_RB01",
                                                              "barcode02_RB02"],
                                                             ["podning", "Væv"],
                                                             ["abundance",
                                                              "abundance"]],
                                                            names = ["prøvenummer",
                                                                     "prøvemateriale", None])
        expected_abundance = pd.DataFrame(data = expected_abundance_values, index = expected_index,
                                          columns = expected_abundance_cols)
        expected_count_values = [[4, np.nan],
                                 [15, 16],
                                 [np.nan, 4],
                                 [1, 0]]
        expected_count_cols = pd.MultiIndex.from_arrays([["barcode01_RB01",
                                                          "barcode02_RB02"],
                                                         ["podning", "Væv"],
                                                         ["counts",
                                                          "counts"]],
                                                        names = ["prøvenummer",
                                                                 "prøvemateriale", None])
        expected_count = pd.DataFrame(data = expected_count_values, index = expected_index,
                                      columns = expected_count_cols)
        test_sheet = (pathlib.Path(__file__).parent / "data" / "summarize_emu"
                      / "test_results_mads.xlsx")
        summarize_emu.write_to_sheets(expected_merged, test_sheet)
        test_merged = pd.read_excel(test_sheet, sheet_name = "overview", index_col = 0,
                                    header = [0, 1, 2])
        test_merged.loc[["Placeholderia bielefeldensis",
                         "Placeholderia fakeorum",
                         "unassigned"], pd.IndexSlice[["barcode01_RB01"], :,
                                        ["med"]]] = ""
        test_merged.loc[["Placeholderia fakeorum",
                         "Placeholderia testfacei",
                         "unassigned"], pd.IndexSlice[["barcode02_RB02"], :,
                                        ["med"]]] = ""
        test_abundance = pd.read_excel(test_sheet, sheet_name = "abundance", index_col = 0,
                                       header = [0, 1, 2])
        test_count = pd.read_excel(test_sheet, sheet_name = "count", index_col = 0,
                                   header = [0, 1, 2])
        pd.testing.assert_frame_equal(test_merged, expected_merged, check_dtype = False)
        pd.testing.assert_frame_equal(test_abundance, expected_abundance, check_dtype = False)
        pd.testing.assert_frame_equal(test_count, expected_count, check_dtype = False)

    def test_handle_blank_anatomy(self):
        """Handle a blank non-sample field in the multiindex."""
        expected_values = [[20.00, 4, "", np.nan, np.nan, np.nan],
                           [75.00, 15, "", 80.00, 16, ""],
                           [np.nan, np.nan, np.nan, 20.00, 4, ""],
                           [5.00, 1, "", 0.00, 0, ""]]
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
                                                      ["podning",
                                                       "podning",
                                                       "podning",
                                                       "",
                                                       "",
                                                       ""],
                                                      ["abundance", "counts",
                                                       "med",
                                                       "abundance", "counts",
                                                       "med"]],
                                                     names = ["prøvenummer",
                                                              "prøvemateriale", None])
        expected_merged = pd.DataFrame(data = expected_values, index = expected_index,
                                       columns = expected_columns)
        expected_abundance_values = [[20.00, np.nan],
                                     [75.00, 80.00],
                                     [np.nan, 20.00],
                                     [5.00, 0.00]]
        expected_abundance_cols = pd.MultiIndex.from_arrays([["barcode01_RB01",
                                                              "barcode02_RB02"],
                                                             ["podning", ""],
                                                             ["abundance",
                                                              "abundance"]],
                                                            names = ["prøvenummer",
                                                                     "prøvemateriale", None])
        expected_abundance = pd.DataFrame(data = expected_abundance_values, index = expected_index,
                                          columns = expected_abundance_cols)
        expected_count_values = [[4, np.nan],
                                 [15, 16],
                                 [np.nan, 4],
                                 [1, 0]]
        expected_count_cols = pd.MultiIndex.from_arrays([["barcode01_RB01",
                                                          "barcode02_RB02"],
                                                         ["podning", ""],
                                                         ["counts",
                                                          "counts"]],
                                                        names = ["prøvenummer",
                                                                 "prøvemateriale", None])
        expected_count = pd.DataFrame(data = expected_count_values, index = expected_index,
                                      columns = expected_count_cols)
        test_sheet = (pathlib.Path(__file__).parent / "data" / "summarize_emu"
                      / "test_results_mads_blank.xlsx")
        summarize_emu.write_to_sheets(expected_merged, test_sheet)
        test_merged = pd.read_excel(test_sheet, sheet_name = "overview", index_col = 0,
                                    header = [0, 1, 2])
        # if the header was blank then it'll be renamed to "Unnamed [n]" - handle all of these
        test_merged = test_merged.rename(columns = lambda colname: "" if "Unnamed" in colname
                                                                   else colname)
        test_merged.loc[["Placeholderia bielefeldensis",
                         "Placeholderia fakeorum",
                         "unassigned"], pd.IndexSlice[["barcode01_RB01"], :,
                                        ["med"]]] = ""
        test_merged.loc[["Placeholderia fakeorum",
                         "Placeholderia testfacei",
                         "unassigned"], pd.IndexSlice[["barcode02_RB02"], :,
                                        ["med"]]] = ""
        test_abundance = pd.read_excel(test_sheet, sheet_name = "abundance", index_col = 0,
                                       header = [0, 1, 2])
        test_abundance = test_abundance.rename(columns = lambda colname: "" if "Unnamed" in colname
                                                                         else colname)
        test_count = pd.read_excel(test_sheet, sheet_name = "count", index_col = 0,
                                   header = [0, 1, 2])
        test_count = test_count.rename(columns = lambda colname: "" if "Unnamed" in colname
                                                                 else colname)
        pd.testing.assert_frame_equal(test_merged, expected_merged)
        pd.testing.assert_frame_equal(test_abundance, expected_abundance)
        pd.testing.assert_frame_equal(test_count, expected_count)


class TestGetNanostatCount(unittest.TestCase):
    def test_fail_file_not_found(self):
        """Fail if the Nanostat result file does not exist."""
        nanostat_file = pathlib.Path(__file__).parent / "data" / "summarize_emu" / "no_such_file"
        error_msg = f"Nanostat report {nanostat_file} not found."
        with pytest.raises(FileNotFoundError, match=re.escape(error_msg)):
            summarize_emu.extract_nanostat_read_count(nanostat_file)

    def test_fail_bad_format(self):
        """Fail informatively for an invalid file format."""
        nanostat_file = (pathlib.Path(__file__).parent / "data" / "summarize_emu"
                         / "nanostat_bad_format.tsv")
        error_msg = "Nanostat report is malformed - cannot find number_of_reads."
        with pytest.raises(KeyError, match = re.escape(error_msg)):
            summarize_emu.extract_nanostat_read_count(nanostat_file)

    def test_success_get_read_count(self):
        """Successfully extract the read count from the Nanostat result file."""
        nanostat_file = (pathlib.Path(__file__).parent / "data" / "summarize_emu"
                         / "nanostat_success.tsv")
        expected_reads = 6176
        test_reads = summarize_emu.extract_nanostat_read_count(nanostat_file)
        assert test_reads == expected_reads


class TestGetKrakenReadCount(unittest.TestCase):
    def test_fail_file_not_found(self):
        """Fail if the Kraken result file does not exist."""
        kraken_file = pathlib.Path(__file__).parent / "data" / "summarize_emu" / "no_such_file"
        error_msg = f"Kraken report {kraken_file} not found."
        with pytest.raises(FileNotFoundError, match = re.escape(error_msg)):
            summarize_emu.get_kraken_read_stats(kraken_file)

    def test_fail_bad_format(self):
        """Fail informatively for an invalid file format."""
        kraken_file = (pathlib.Path(__file__).parent / "data" / "summarize_emu"
                         / "kraken_bad_format.tsv")
        error_msg = "Kraken report is malformed. Check that you are supplying a --report file."
        with pytest.raises(KeyError, match = re.escape(error_msg)):
            summarize_emu.get_kraken_read_stats(kraken_file)

    def test_fail_other_error(self):
        """Don't override other value errors."""
        kraken_file = (pathlib.Path(__file__).parent / "data" / "summarize_emu"
                     / "kraken_bad_encoding.tsv")  # IBM500 encoded to make utf-8 reader complain
        error_msg = "'utf-8' codec can't decode"
        with pytest.raises(ValueError, match = re.escape(error_msg)):
            summarize_emu.get_kraken_read_stats(kraken_file)

    def test_success_get_read_count(self):
        """Successfully extract the read counts from the Kraken result file."""
        kraken_file = (pathlib.Path(__file__).parent / "data" / "summarize_emu"
                       / "kraken_success.tsv")
        expected_results = pd.DataFrame(data={"human": [10], "remaining": [90], "total": [100]})
        test_results = summarize_emu.get_kraken_read_stats(kraken_file)
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_success_log_no_human(self):
        """Log if no human reads were removed."""
        kraken_file = (pathlib.Path(__file__).parent / "data" / "summarize_emu"
                       / "kraken_no_human.tsv")
        expected_results = pd.DataFrame(data = {"human": [0], "remaining": [100], "total": [100]})
        log_msg = f"INFO:summarize_emu:No human reads reported in {kraken_file}."
        with self.assertLogs("summarize_emu", level="INFO") as logged:
            test_results = summarize_emu.get_kraken_read_stats(kraken_file)
            assert log_msg in logged.output
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_success_log_no_remaining(self):
        """Log if no reads remain after human depletion."""
        kraken_file = (pathlib.Path(__file__).parent / "data" / "summarize_emu"
                       / "kraken_no_remaining.tsv")
        expected_results = pd.DataFrame(data = {"human": [100], "remaining": [0], "total": [100]})
        log_msg = f"INFO:summarize_emu:No remaining reads reported in {kraken_file}."
        with self.assertLogs("summarize_emu", level = "INFO") as logged:
            test_results = summarize_emu.get_kraken_read_stats(kraken_file)
            assert log_msg in logged.output
        pd.testing.assert_frame_equal(expected_results, test_results)

    def test_success_empty(self):
        """Log if the file contains no reads."""
        kraken_file = (pathlib.Path(__file__).parent / "data" / "summarize_emu"
                       / "kraken_empty.tsv")
        expected_results = pd.DataFrame(data = {"human": [0], "remaining": [0], "total": [0]})
        no_human_msg = f"INFO:summarize_emu:No human reads reported in {kraken_file}."

        no_bact_msg = f"INFO:summarize_emu:No remaining reads reported in {kraken_file}."
        with self.assertLogs("summarize_emu", level = "INFO") as logged:
            test_results = summarize_emu.get_kraken_read_stats(kraken_file)
            assert no_human_msg in logged.output
            assert no_bact_msg in logged.output
        pd.testing.assert_frame_equal(expected_results, test_results)