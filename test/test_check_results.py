import io
import logging
import pathlib
import re
import unittest

from unittest import mock
from unittest.mock import PropertyMock

import pandas as pd
import pytest

import check_results
import version


class TestCheckFilePresence(unittest.TestCase):
    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog
    
    def test_success(self):
        """Successfully find all relevant files."""
        test_dir = pathlib.Path(__file__).parent / "data"/"check_results"/"success_emu_dir"
        check_files = check_results.check_files_present(test_dir)
        assert check_files

    def test_missing_emu_files(self):
        """Alert when the Emu report is missing"""
        test_dir = pathlib.Path(__file__).parent / "data" / "check_results" / "empty_dir"
        warn_msg = "Emu report is missing. Cannot evaluate Emu results."
        with self._caplog.at_level(logging.WARNING, logger = "QATest"):
            check_files = check_results.check_files_present(test_dir)
            assert ("QATest", logging.WARNING, warn_msg) in self._caplog.record_tuples
        assert not check_files

    def test_missing_raw_backup(self):
        """Alert when the backup file is missing."""
        test_dir = pathlib.Path(__file__).parent / "data" / "check_results" / "emu_dir_no_backup"
        warn_msg = "WARNING:QATest:Raw TSV backup of Emu report is missing."
        with self.assertLogs("QATest") as logged:
            check_files = check_results.check_files_present(test_dir)
            assert warn_msg in logged.output
        assert not check_files

    def test_too_many_emus(self):
        """Alert when there are multiple Emu reports"""
        test_dir = pathlib.Path(__file__).parent / "data" / "check_results" / "multiple_emus"
        warn_msg = ("WARNING:QATest:Multiple Emu summaries found, need only one. "
                    "Cannot evaluate Emu results.")
        with self.assertLogs("QATest") as logged:
            check_files = check_results.check_files_present(test_dir)
            assert warn_msg in logged.output
        assert not check_files

@mock.patch(f"{check_results.__name__}.__version__")
class TestCheckEmuResults(unittest.TestCase):
    def test_fail_version_mismatch(self, mock_version):
        """Fail if a report is being checked with a different version of the pipeline than the
        one that generated it."""
        mock_version.__str__.return_value = "0.4.2"
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_bad_version_emu-combined.xlsx")
        error_msg = ("Report was created with Version_1.2.3, is being checked with Version_0.4.2.\n"
                     "Cannot check reports from a different version of the pipeline.")
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            check_results.check_emu_result_file(test_report)

    def test_all_good(self, mock_version):
        """Don't complain for an Emu report without any issues."""
        mock_version.__str__.return_value = "0.4.2"
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results" / "success_emu_dir"
                       / "RUN0001_emu-combined.xlsx")
        check_report = check_results.check_emu_result_file(test_report)
        assert check_report

    def test_success_minor_abundance_diff(self, mock_version):
        """Report success if there is a small difference in abundance in the positive control
        (<0.5 percent points)."""
        mock_version.__str__.return_value = "0.4.2"
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_minor_diff_emu-combined.xlsx")
        check_report = check_results.check_emu_result_file(test_report)
        assert check_report

    def test_warn_missing_tabs(self, mock_version):
        """Complain if one or more of the output sheets are missing."""
        mock_version.__str__.return_value = "0.4.2"
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_missing_tab_emu-combined.xlsx")
        warn_msg = ("WARNING:QATest:Tab(s) ['abundance', 'count'] not found in Emu report. "
                    "Cannot evaluate results for these tabs.")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_emu_result_file(test_report)
            assert warn_msg in logged.output
        assert not check_report

    def test_warn_wrong_header(self, mock_version):
        """Warn if the header rows don't match what is expected."""
        mock_version.__str__.return_value = "0.4.2"
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_bad_header_emu-combined.xlsx")
        mismatch_header = pd.MultiIndex.from_arrays([["anatomi", "anatomi"],
                                                     ["expected", "found"]])
        mismatch_index = pd.Index(["F99123456"], name="prøvenr")
        mismatched = pd.DataFrame(data=[["Svælg/tonsil", "Næse"]], index = mismatch_index,
                                  columns = mismatch_header)
        warn_msg = ("WARNING:QATest:Sample metadata differ from expected sample metadata in tab"
                    " overview:\n"
                    f"{mismatched.to_string()}")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_emu_result_file(test_report)
            assert warn_msg in logged.output
        assert not check_report

    def test_warn_broken_header(self, mock_version):
        """Warn if the format of the header differs from the expected format."""
        mock_version.__str__.return_value = "0.4.2"
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_bad_cols_emu-combined.xlsx")
        cols_per_sample = 3
        expected_index = pd.Index([*["NegK_Sanger"] * cols_per_sample,
                                   *["PosK"] * cols_per_sample,
                                   *["F99123457"] * cols_per_sample,
                                   *["F99123456"] * cols_per_sample,
                                   *["F99123458"] * cols_per_sample
                                   ],
                                  name = "prøvenr")
        expected_cols = pd.Index(["run",
                                  "pipeline_version",
                                  "barcode",
                                  "modtagedato",
                                  "patient",
                                  "prøvemateriale",
                                  "anatomi",
                                  "indikation",
                                  "total_before_qc",
                                  "total_after_qc",
                                  "human"])
        found_index = pd.Index([*["PosK"] * cols_per_sample,
                                *["NegK_Sanger"] * cols_per_sample,
                                *["F99123457"] * cols_per_sample,
                                *["F99123456"] * cols_per_sample,
                                *["F99123458"] * cols_per_sample
                                ],
                               name = "prøvenr")
        found_cols = pd.Index(["run",
                               "pipeline_version",
                               "barcode",
                               "modtagedato",
                               "patient",
                               "prøvemateriale",
                               "anatomi",
                               "indikation",
                               "total_before_qc",
                               "total_after_qc",
                               "human"])
        warn_msg = ("WARNING:QATest:Sample metadata labels differ from expected sample metadata"
                    " - could not compare. \n"
                    f"Expected index: {expected_index}\n"
                    f"Found index: {found_index}\n"
                    f"Expected columns: {expected_cols}\n"
                    f"Found columns: {found_cols}")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_emu_result_file(test_report)
            print(logged.output)
            assert warn_msg in logged.output
        assert not check_report

    def test_warn_wrong_organism_main_tab(self, mock_version):
        """Warn if one of the test samples isn't the organism we expect it to be (overview tab)."""
        mock_version.__str__.return_value = "0.4.2"
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_bad_organism_main_emu-combined.xlsx")
        mismatch_header = pd.MultiIndex.from_arrays([["organism", "organism"],
                                                     ["expected", "found"]])
        expected_data = pd.DataFrame(data=[["Streptococcus agalactiae",
                                            "Placeholderia bielefeldensis"]],
                                     index = pd.Index(["F99123457"], name="prøvenr"),
                                     columns=mismatch_header)
        print(expected_data)
        print(expected_data.to_string())
        warn_msg = ("WARNING:QATest:Incorrect organism for one or more samples."
                    " Expected organism(s):\n"
                    f"{expected_data.to_string()}")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_emu_result_file(test_report)
            print(warn_msg)
            print(logged.output)
            assert warn_msg in logged.output
        assert not check_report

    def test_warn_wrong_positive_control_main(self, mock_version):
        """Warn if the positive control does not contain the expected species (overview tab)."""
        mock_version.__str__.return_value = "0.4.2"
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_bad_positive_control_main_emu-combined.xlsx")
        warn_msg = ("WARNING:QATest:Positive control should contain "
                    "['Bacillus subtilis', "
                    "'Enterococcus faecalis', "
                    "'Escherichia coli', "
                    "'Limosilactobacillus fermentum', "
                    "'Listeria monocytogenes', "
                    "'Pseudomonas aeruginosa', "
                    "'Salmonella enterica', "
                    "'Staphylococcus aureus'], "
                    "contains ['Bacillus subtilis', "
                    "'Enterococcus faecalis', "
                    "'Escherichia coli', "
                    "'Limosilactobacillus fermentum', "
                    "'Listeria monocytogenes', "
                    "'Placeholderia bielefeldensis', "
                    "'Salmonella enterica', "
                    "'Staphylococcus aureus']"
                    " (missing: {'Pseudomonas aeruginosa'}, "
                    "extra: {'Placeholderia bielefeldensis'}")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_emu_result_file(test_report)
            assert warn_msg in logged.output
        assert not check_report

    def test_warn_wrong_abundance_main(self, mock_version):
        """Warn if the relative abundance of positive control organisms varies too much
        from what is expected (overview tab)."""
        mock_version.__str__.return_value = "0.4.2"
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_bad_abundance_main_emu-combined.xlsx")
        expected_data = pd.DataFrame(data = {"expected": [19.17],
                                             "found": [25.00]},
                                     index = pd.Index(["Salmonella enterica"], name = "organism"))
        warn_msg = ("WARNING:QATest:Different abundance in positive control for "
                    "['Salmonella enterica']."
                    " Expected abundance:\n"
                    f"{expected_data.to_string()}")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_emu_result_file(test_report)
            assert warn_msg in logged.output
        assert not check_report

    def test_multiple_mismatches_one_sample(self, mock_version):
        """Report on multiple issues with one samples."""
        mock_version.__str__.return_value = "0.4.2"
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_bad_positive_control_abundance_main_emu-combined.xlsx")
        wrong_organism = ("WARNING:QATest:Positive control should contain "
                          "['Bacillus subtilis', "
                          "'Enterococcus faecalis', "
                          "'Escherichia coli', "
                          "'Limosilactobacillus fermentum', "
                          "'Listeria monocytogenes', "
                          "'Pseudomonas aeruginosa', "
                          "'Salmonella enterica', "
                          "'Staphylococcus aureus'], "
                          "contains ['Bacillus subtilis', "
                          "'Enterococcus faecalis', "
                          "'Escherichia coli', "
                          "'Limosilactobacillus fermentum', "
                          "'Listeria monocytogenes', "
                          "'Placeholderia bielefeldensis', "
                          "'Salmonella enterica', "
                          "'Staphylococcus aureus']"
                          " (missing: {'Pseudomonas aeruginosa'}, "
                          "extra: {'Placeholderia bielefeldensis'}")
        expected_data = pd.DataFrame(data = {"expected": [19.17],
                                             "found": [25.00]},
                                     index = pd.Index(["Salmonella enterica"], name = "organism"))
        wrong_abundance = ("WARNING:QATest:Different abundance in positive control for "
                           "['Salmonella enterica']."
                           " Expected abundance:\n"
                           f"{expected_data.to_string()}")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_emu_result_file(test_report)
            assert wrong_organism in logged.output
            assert wrong_abundance in logged.output
        assert not check_report

    def test_multi_sample_mismatches(self, mock_version):
        """Report on multiple issues with multiple samples."""
        mock_version.__str__.return_value = "0.4.2"
        test_report = (pathlib.Path(__file__).parent / "data" / "check_results"
                       / "RUN0001_multi_fail_main_emu-combined.xlsx")
        mismatch_header = pd.MultiIndex.from_arrays([["organism", "organism"],
                                                     ["expected", "found"]])
        expected_data = pd.DataFrame(data = [["Cutibacterium acnes",
                                              "Placeholderia bielefeldensis"],
                                             ["Streptococcus agalactiae",
                                              "Placeholderia bielefeldensis"]],
                                     index = pd.Index(["F99123456", "F99123457"],
                                                      name = "prøvenr"),
                                     columns = mismatch_header)
        warn_msg = ("WARNING:QATest:Incorrect organism for one or more samples."
                    " Expected organism(s):\n"
                    f"{expected_data.to_string()}")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_emu_result_file(test_report)
            assert warn_msg in logged.output
        assert not check_report


@mock.patch(f"{check_results.__name__}.__version__")
class TestCheckAllQC(unittest.TestCase):
    def test_warn_missing_files(self, mock_version):
        """Warn if there are missing files."""
        mock_version.__str__.return_value = "0.4.2"
        test_dir = pathlib.Path(__file__).parent / "data" / "check_results" / "empty_dir"
        warn_msg = "WARNING:QATest:Emu report is missing. Cannot evaluate Emu results."
        with self.assertLogs("QATest") as logged:
            check_files = check_results.check_all_qc(test_dir)
            assert warn_msg in logged.output
        assert not check_files

    def test_warn_wrong_results(self, mock_version):
        """Warn if any results differ from the expected results."""
        mock_version.__str__.return_value = "0.4.2"
        test_dir = pathlib.Path(__file__).parent / "data" / "check_results" / "bad_emu_dir"
        expected_data = pd.DataFrame(data = {"expected": [19.17],
                                             "found": [25.00]},
                                     index = pd.Index(["Salmonella enterica"], name = "organism"))
        warn_msg = ("WARNING:QATest:Different abundance in positive control for "
                    "['Salmonella enterica']."
                    " Expected abundance:\n"
                    f"{expected_data.to_string()}")
        with self.assertLogs("QATest") as logged:
            check_report = check_results.check_all_qc(test_dir)
            assert warn_msg in logged.output
        assert not check_report

    def test_success(self, mock_version):
        """Report success if all checks pass."""
        mock_version.__str__.return_value = "0.4.2"
        test_dir = pathlib.Path(__file__).parent / "data"/"check_results"/"success_emu_dir"
        success_msg = "INFO:QATest:All QC checks passed"
        with self.assertLogs("QATest") as logged:
            check_files = check_results.check_all_qc(test_dir)
            assert success_msg in logged.output
        assert check_files

@mock.patch(f"{check_results.__name__}.__version__")
class TestCheckResults(unittest.TestCase):
    log_format = re.compile(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3}"
                            r" - QATest - (INFO|WARNING|ERROR) - [-:a-zA-Z0-9_./# ]+")
    def tearDown(self):
        """Clean up existing logfiles."""
        (pathlib.Path(__file__).parent / "data" / "check_results"
         / "empty_dir" / "logs" / "pipeline_qa.log").unlink(missing_ok = True)
        (pathlib.Path(__file__).parent / "data" / "check_results"
         / "success_emu_dir" / "logs" / "pipeline_qa.log").unlink(missing_ok = True)
        (pathlib.Path(__file__).parent / "data" / "check_results"
         / "custom.log").unlink(missing_ok = True)

    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    def test_warn_missing_results(self, mock_version):
        """Warn if any results are missing, and log this properly."""
        mock_version.__str__.return_value = "0.4.2"
        test_dir = pathlib.Path(__file__).parent / "data" / "check_results" / "empty_dir"
        expected_logfile = test_dir / "logs" / "pipeline_qa.log"
        error_msg = "One or more QC steps failed. Check log for details."
        warn_msg = "Emu report is missing. Cannot evaluate Emu results."
        assert not expected_logfile.exists()
        with pytest.raises(SystemExit, match="1"), self._caplog.at_level(logging.INFO,
                                                                         logger="QATest"):
            check_results.check_results([str(test_dir)])
        assert ("QATest", logging.ERROR, error_msg) in self._caplog.record_tuples
        assert ("QATest", logging.WARNING, warn_msg) in self._caplog.record_tuples
        assert expected_logfile.exists()
        with open(expected_logfile, "r") as read_log:
            log_data = read_log.readlines()
        for line in log_data:
            assert re.match(self.log_format, line.strip())


    def test_success(self, mock_version):
        """Report success if all checks pass, and log this properly."""
        mock_version.__str__.return_value = "0.4.2"
        test_dir = pathlib.Path(__file__).parent / "data" / "check_results" / "success_emu_dir"
        expected_logfile = (pathlib.Path(__file__).parent / "data" / "check_results"
         / "custom.log")
        success_msg = "All QC checks passed"
        assert not expected_logfile.exists()
        with pytest.raises(SystemExit, match = "0"), self._caplog.at_level(logging.INFO,
                                                                           logger = "QATest"):
            check_results.check_results([str(test_dir), "--logfile", str(expected_logfile)])
        assert ("QATest", logging.INFO, success_msg) in self._caplog.record_tuples
        assert expected_logfile.exists()
        with open(expected_logfile, "r") as read_log:
            log_data = read_log.readlines()
        for line in log_data:
            assert re.match(self.log_format, line.strip())


    def test_different_logfile(self, mock_version):
        """Set different logfile."""
        mock_version.__str__.return_value = "0.4.2"
        test_dir = pathlib.Path(__file__).parent / "data" / "check_results" / "success_emu_dir"
        expected_logfile = test_dir / "logs" / "pipeline_qa.log"
        success_msg = "All QC checks passed"
        assert not expected_logfile.exists()
        with pytest.raises(SystemExit, match = "0"), self._caplog.at_level(logging.INFO,
                                                                           logger = "QATest"):
            check_results.check_results([str(test_dir)])
        assert ("QATest", logging.INFO, success_msg) in self._caplog.record_tuples
        assert expected_logfile.exists()
        with open(expected_logfile, "r") as read_log:
            log_data = read_log.readlines()
        for line in log_data:
            assert re.match(self.log_format, line.strip())