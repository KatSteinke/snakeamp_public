import pathlib
import unittest

import pandas as pd
import pytest

import check_results


class TestCheckFilePresence(unittest.TestCase):
    def test_success(self):
        """Successfully find all relevant files."""
        test_dir = pathlib.Path(__file__).parent / "data"/"check_results"/"success_emu_dir"
        check_files = check_results.check_files_present(test_dir)
        assert check_files

    def test_missing_emu_files(self):
        """Alert when the Emu report is missing"""
        test_dir = pathlib.Path(__file__).parent / "data" / "check_results" / "empty_dir"
        warn_msg = "WARNING:QATest:Emu report is missing. Cannot evaluate Emu results."
        with self.assertLogs("QATest") as logged:
            check_files = check_results.check_files_present(test_dir)
            assert warn_msg in logged.output
        assert not check_files


class TestCheckEmuResults(unittest.TestCase):
    def test_all_good(self):
        """Don't complain for an Emu report without any issues."""
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results" / "success_emu_dir"
                       / "RUN0001_emu-combined.xlsx")
        check_report = check_results.check_emu_result_file(test_report)
        assert check_report

    def test_warn_missing_tabs(self):
        """Complain if one or more of the output sheets are missing."""
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_missing_tab_emu-combined.xlsx")
        warn_msg = ("WARNING:QATest:Tab(s) ['abundance', 'count'] not found in Emu report. "
                    "Cannot evaluate results for these tabs.")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_emu_result_file(test_report)
            assert warn_msg in logged.output
        assert not check_report

    def test_warn_wrong_header(self):
        """Warn if the header rows don't match what is expected."""
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_bad_header_emu-combined.xlsx")
        mismatch_header = pd.MultiIndex.from_arrays([["anatomi", "anatomi"],
                                                     "expected", "found"])
        mismatch_index = pd.Index(["1199123456-1"], name="prøvenr")
        mismatched = pd.DataFrame(data=[["Svælg/tonsil", "Næse"]], index = mismatch_index,
                                  columns = mismatch_header)
        warn_msg = ("WARNING:QATest:Sample metadata differ from expected sample metadata in tab"
                    " 'overview':\n"
                    f"{mismatched.to_string()}")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_emu_result_file(test_report)
            assert warn_msg in logged.output
        assert not check_report

    def test_warn_wrong_organism_main_tab(self):
        """Warn if one of the test samples isn't the organism we expect it to be (overview tab)."""
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_bad_organism_main_emu-combined.xlsx")
        expected_data = pd.DataFrame(data={"expected": ["Propionibacterium acnes"],
                                           "found": ["Placeholderia bielefeldensis"]},
                                     index = pd.Index(["1199123456-1"], name="prøvenr"))
        warn_msg = ("WARNING:QATest:Incorrect organism for sample number(s) ['1199123456-1']."
                    " Expected organism(s):\n"
                    f"{expected_data.to_string()}")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_emu_result_file(test_report)
            assert warn_msg in logged.output
        assert not check_report

    def test_warn_wrong_positive_control_main(self):
        """Warn if the positive control does not contain the expected species (overview tab)."""
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_bad_positive_control_main_emu-combined.xlsx")
        warn_msg = ("WARNING:QATest:Positive control should contain ['Bacillus subtilis',"
                    " 'Staphylococcus aureus', 'Listeria monocytogenes', 'Salmonella enterica',"
                    " 'Escherichia coli', 'Enterococcus faecalis', 'Pseudomonas aeruginosa'],"
                    "contains ['Bacillus subtilis',"
                    " 'Staphylococcus aureus', 'Listeria monocytogenes', 'Salmonella enterica',"
                    " 'Escherichia coli', 'Enterococcus faecalis', 'Placeholderia bielefeldensis']"
                    "(missing: {'Pseudomonas aeruginosa'}, extra: {'Placeholderia bielefeldensis'}")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_emu_result_file(test_report)
            assert warn_msg in logged.output
        assert not check_report

    def test_warn_wrong_abundance_main(self):
        """Warn if the relative abundance of positive control organisms varies too much
        from what is expected (overview tab)."""
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_bad_abundance_main_emu-combined.xlsx")
        expected_data = pd.DataFrame(data = {"expected": [20.04],
                                             "found": [25]},
                                     index = pd.Index(["Salmonella enterica"], name = "prøvenr"))
        warn_msg = ("WARNING:QATest:Different abundance in positive control for "
                    "['Salmonella enterica']."
                    " Expected abundance:\n"
                    f"{expected_data.to_string()}")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_emu_result_file(test_report)
            assert warn_msg in logged.output
        assert not check_report

    def test_multiple_mismatches_one_sample(self):
        """Report on multiple issues with one samples."""
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_bad_positive_control_abundance_main_emu-combined.xlsx")
        wrong_organism = ("WARNING:QATest:Positive control should contain ['Bacillus subtilis',"
                    " 'Staphylococcus aureus', 'Listeria monocytogenes', 'Salmonella enterica',"
                    " 'Escherichia coli', 'Enterococcus faecalis', 'Pseudomonas aeruginosa'],"
                    "contains ['Bacillus subtilis',"
                    " 'Staphylococcus aureus', 'Listeria monocytogenes', 'Salmonella enterica',"
                    " 'Escherichia coli', 'Enterococcus faecalis', 'Placeholderia bielefeldensis']"
                    "(missing: {'Pseudomonas aeruginosa'}, extra: {'Placeholderia bielefeldensis'}")
        expected_data = pd.DataFrame(data = {"expected": [20.04],
                                             "found": [25]},
                                     index = pd.Index(["Salmonella enterica"], name = "prøvenr"))
        wrong_abundance = ("WARNING:QATest:Different abundance in positive control for "
                    "['Salmonella enterica']."
                    " Expected abundance:\n"
                    f"{expected_data.to_string()}")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_emu_result_file(test_report)
            assert wrong_organism in logged.output
            assert wrong_abundance in logged.output
        assert not check_report

    def test_multi_sample_mismatches(self):
        """Report on multiple issues with multiple samples."""
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_bad_organism_main_emu-combined.xlsx")
        expected_data = pd.DataFrame(data = {"expected": ["Propionibacterium acnes",
                                                          "Streptococcus agalactiae"],
                                             "found": ["Placeholderia bielefeldensis",
                                                       "Placeholderia bielefeldensis"]},
                                     index = pd.Index(["1199123456-1", "1199123457-1"],
                                                      name = "prøvenr"))
        warn_msg = ("WARNING:QATest:Incorrect organism for sample number(s) ['1199123456-1',"
                    " '1199123457-1']."
                    " Expected organism(s):\n"
                    f"{expected_data.to_string()}")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_emu_result_file(test_report)
            assert warn_msg in logged.output
        assert not check_report

