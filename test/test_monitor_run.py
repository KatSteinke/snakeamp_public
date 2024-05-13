import pathlib
import re
import unittest

from unittest import mock

import pytest

import monitor_run


class TestWaitForFile(unittest.TestCase):
    active_config = {"sample_number_settings": {"sample_number_format":
                                                    '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})-\d?',
                                                "sample_numbers_in": "letter",
                                                "sample_numbers_out": "letter",
                                                "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                "number_to_letter": {"70": "P",
                                                                     "30": "B",
                                                                     "10": "D",
                                                                     "50": "T"},
                                                "date_settings":
                                                    {"splice_in_date": False,
                                                     "length_without_date": 8,
                                                     "splice_after": 2},
                                                "negative_control": 'NegK[a-zA-Z0-9]*',
                                                "positive_control": {}},
                     "barcode_format": "NB[0-9]{2}",  # format of barcodes in runsheet
                     "barcode_prefix": "NB",
                     "amplicon_type": "16S",
                     "debug": False}
    config_path = pathlib.Path(__file__).parent / "data"/"monitor_run"/"test_config.yaml"
    test_run = monitor_run.AmpliconRun(sequence_dir = pathlib.Path(__file__).parent / "data"
                                                      /"monitor_run"/"miniondir"/"test1",
                                       outdir = pathlib.Path(__file__).parent / "data"
                                                / "monitor_run" / "test_outdir",
                                       runsheet = pathlib.Path(__file__).parent / "data"
                                                  / "utilities_test"
                                                  / "test_nanopore_runsheet_16s_only.xlsx",
                                       active_config = active_config,
                                       configfile = config_path,
                                       test_run = False)

    def test_timeout(self):
        """Stop monitoring when the timeout has been reached."""
        error_msg = ("No file matching pattern test_summary*.txt "
                     f"found in {self.test_run.sequence_dir.parent}")
        with pytest.raises(FileNotFoundError, match=re.escape(error_msg)):
            monitor_run.start_on_file_found(self.test_run, "test_summary*.txt", watch_interval = 1,
                                            watch_timeout = 5)

    def test_fail_parent_dir_not_found(self):
        """Don't start monitoring if the parent directory does not exist."""
        test_run = monitor_run.AmpliconRun(sequence_dir = pathlib.Path(__file__).parent / "data"
                                                          / "monitor_run" / "miniondir" / "test2",
                                           outdir = pathlib.Path(__file__).parent / "data"
                                                    / "monitor_run" / "test_outdir",
                                           runsheet = pathlib.Path(__file__).parent / "data"
                                                      / "utilities_test"
                                                      / "test_nanopore_runsheet_16s_only.xlsx",
                                           active_config = self.active_config,
                                           configfile = self.config_path,
                                           test_run = False)
        error_msg = f"Parent directory {test_run.sequence_dir} does not exist."
        with pytest.raises(FileNotFoundError, match = re.escape(error_msg)):
            monitor_run.start_on_file_found(test_run, "*/final_summary*.txt",
                                            watch_interval = 1,
                                            watch_timeout = 5, dry_run = True)

    def test_fail_bad_parent_dir(self):
        """Don't start monitoring if the parent directory does not have the expected structure."""
        test_run = monitor_run.AmpliconRun(sequence_dir = pathlib.Path(__file__).parent / "data"
                                                          / "monitor_run" / "miniondir" / "test3",
                                           outdir = pathlib.Path(__file__).parent / "data"
                                                    / "monitor_run" / "test_outdir",
                                           runsheet = pathlib.Path(__file__).parent / "data"
                                                      / "utilities_test"
                                                      / "test_nanopore_runsheet_16s_only.xlsx",
                                           active_config = self.active_config,
                                           configfile = self.config_path,
                                           test_run = False)
        error_msg = (f"Parent directory {test_run.sequence_dir} "
                     "does not contain a rawdata directory.")
        with pytest.raises(FileNotFoundError, match = re.escape(error_msg)):
            monitor_run.start_on_file_found(test_run, "*/final_summary*.txt",
                                            watch_interval = 1,
                                            watch_timeout = 5, dry_run = True)

    def test_success_dryrun(self):
        """Successfully print the pipeline start command in dry run mode."""
        expected_command = ['echo', (f'"nomad job dispatch -meta indir={self.test_run.sequence_dir}'
                                     f' -meta outdir={self.test_run.outdir} '
                                     f'-meta runsheet={self.test_run.runsheet} '
                                     '16s-snake-emu-prod '
                                     f'{str(self.test_run.configfile)}"')]
        expected_log = ("INFO:launch_run:Found */final_summary*.txt "
                        f"in {self.test_run.sequence_dir / 'rawdata'} after 0 seconds.")
        with self.assertLogs("launch_run", level = "INFO") as logged:
            test_command = monitor_run.start_on_file_found(self.test_run, "*/final_summary*.txt",
                                                           dry_run = True, watch_interval = 10,
                                                           watch_timeout = 30).args
            assert expected_log in logged.output
        assert test_command == expected_command

    # we don't want to run the actual nomad command when testing...
    @mock.patch(f'{monitor_run.__name__}.get_pipeline_command',
                wraps = monitor_run.get_pipeline_command)
    def test_success_real_run(self, mock_command):
        """Successfully run the output of the pipeline command builder."""
        mock_command.return_value = ["echo", "Hello"]
        expected_command = ["echo", "Hello"]
        test_command = monitor_run.start_on_file_found(self.test_run, "*/final_summary*.txt",
                                                       dry_run = False, watch_interval = 10,
                                                       watch_timeout = 30).args
        assert test_command == expected_command


class TestGetNomadCommand(unittest.TestCase):
    active_config = {"sample_number_settings": {"sample_number_format":
                                                    '([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})-\d?',
                                                "sample_numbers_in": "letter",
                                                "sample_numbers_out": "letter",
                                                "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                "number_to_letter": {"70": "P",
                                                                     "30": "B",
                                                                     "10": "D",
                                                                     "50": "T"},
                                                "date_settings":
                                                    {"splice_in_date": False,
                                                     "length_without_date": 8,
                                                     "splice_after": 2},
                                                "negative_control": 'NegK[a-zA-Z0-9]*',
                                                "positive_control": {}},
                     "barcode_format": "NB[0-9]{2}",  # format of barcodes in runsheet
                     "barcode_prefix": "NB",
                     "amplicon_type": "16S",
                     "debug": False}
    config_path = pathlib.Path(__file__).parent / "data"/"monitor_run"/"test_config.yaml"

    def test_get_routine_command(self):
        indir = pathlib.Path("path/to/indir")
        outdir = pathlib.Path("path/to/outdir")
        runsheet = pathlib.Path("path/to/runsheet")
        seq_run = monitor_run.AmpliconRun(sequence_dir = indir,
                                          outdir = outdir,
                                          runsheet = runsheet,
                                          active_config = self.active_config,
                                          configfile = self.config_path)
        expected_command = ["nomad", "job", "dispatch",
                            "-meta", "indir=path/to/indir",
                            "-meta", "outdir=path/to/outdir",
                            "-meta", "runsheet=path/to/runsheet",
                            "16s-snake-emu-prod", str(self.config_path)]
        test_command = monitor_run.get_pipeline_command(seq_run)
        assert expected_command == test_command

    def test_get_test_command(self):
        """Test that the staging version is run if debug is specified (overriding default)."""
        indir = pathlib.Path("path/to/indir")
        outdir = pathlib.Path("path/to/outdir")
        runsheet = pathlib.Path("path/to/runsheet")
        expected_command = ["nomad", "job", "dispatch",
                            "-meta", "indir=path/to/indir",
                            "-meta", "outdir=path/to/outdir",
                            "-meta", "runsheet=path/to/runsheet",
                            "16s-snake-emu-staging", str(self.config_path)]
        seq_run = monitor_run.AmpliconRun(sequence_dir = indir,
                                          outdir = outdir,
                                          runsheet = runsheet,
                                          active_config = self.active_config,
                                          configfile = self.config_path,
                                          test_run = True)
        test_command = monitor_run.get_pipeline_command(seq_run)
        assert expected_command == test_command

    def test_override_test_command(self):
        """Test that a config specifying debug mode can be overridden by the debug flag."""
        indir = pathlib.Path("path/to/indir")
        outdir = pathlib.Path("path/to/outdir")
        runsheet = pathlib.Path("path/to/runsheet")
        test_configfile = pathlib.Path("path/to/config")
        test_config = {"debug": True}
        expected_command = ["nomad", "job", "dispatch",
                            "-meta", "indir=path/to/indir",
                            "-meta", "outdir=path/to/outdir",
                            "-meta", "runsheet=path/to/runsheet",
                            "16s-snake-emu-prod", str(test_configfile)]
        seq_run = monitor_run.AmpliconRun(sequence_dir = indir,
                                          outdir = outdir,
                                          runsheet = runsheet,
                                          active_config = test_config,
                                          configfile = test_configfile,
                                          test_run = False)
        test_command = monitor_run.get_pipeline_command(seq_run)
        assert expected_command == test_command
