import pathlib
import re
import unittest

from collections import namedtuple
from datetime import timedelta
from unittest import mock

import pytest

import monitor_run


class TestCreateAmpliconRun(unittest.TestCase):
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
                     "debug": False,
                     "paths": {"output_base_path": "/path/to/output"},
                     "seq_run_duration_hours": 1}

    def test_create_run(self):
        """Successfully create an AmpliconRun"""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data"/"utilities_test"
                    /"test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        expected_outdir = pathlib.Path("/path/to/output/NANO_Amplicon_Y20990101_RUN0001_XYZ-16S")
        amplicon_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config)
        seq_run_duration = timedelta(hours = float(self.active_config["seq_run_duration_hours"]))
        assert expected_indir == amplicon_run.sequence_dir
        assert runsheet == amplicon_run.runsheet
        assert configfile == amplicon_run.configfile
        assert self.active_config == amplicon_run.active_config
        assert seq_run_duration == amplicon_run.sequencing_time
        assert expected_outdir == amplicon_run.outdir
        assert not amplicon_run.test_run

    def test_initialize_with_seq_time(self):
        """Set the sequencing time when initializing the run."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        new_seq_time = 16
        new_seq_duration = timedelta(hours=new_seq_time)
        amplicon_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config,
                                               sequencing_time = new_seq_time)
        assert amplicon_run.sequencing_time == new_seq_duration

    def test_fail_bad_seq_time(self):
        """Fail when the sequencing time is invalid."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        new_seq_time = -1
        error_msg = "Expected sequencing time must be greater than 0 hours."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            monitor_run.AmpliconRun(sequence_dir = expected_indir,
                                    runsheet = runsheet,
                                    configfile = configfile,
                                    active_config = self.active_config,
                                    sequencing_time = new_seq_time)

    def test_override_test_run(self):
        """Initialize the run as a test run."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        amplicon_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config, test_run = True)
        assert amplicon_run.test_run

    def test_initialize_with_outdir(self):
        """Set a new output directory when initializing the run."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        new_outdir = pathlib.Path("/path/to/outdir")
        amplicon_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config,
                                               outdir = new_outdir)
        assert new_outdir == amplicon_run.outdir

    def test_set_new_outdir(self):
        """Set a new output directory after the run has been initialized."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        initial_outdir = pathlib.Path("/path/to/output/NANO_Amplicon_Y20990101_RUN0001_XYZ-16S")
        amplicon_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config)
        assert initial_outdir == amplicon_run.outdir
        new_outdir = pathlib.Path("/path/to/outdir")
        amplicon_run.outdir = new_outdir
        assert new_outdir == amplicon_run.outdir

    def test_print_run(self):
        """Print the run's attributes in the correct format."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        expected_outdir = pathlib.Path("/path/to/output/NANO_Amplicon_Y20990101_RUN0001_XYZ-16S")
        amplicon_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config)
        expected_repr = (f"AmpliconRun(sequence_dir={expected_indir}, runsheet={runsheet}, "
                         f"configfile={configfile}, active_config={self.active_config}, "
                         f"outdir={expected_outdir},"
                         " sequencing_time="
                         f"{timedelta(hours=self.active_config['seq_run_duration_hours'])})")
        amplicon_repr = amplicon_run.__repr__()
        assert amplicon_repr == expected_repr

    def test_compare_equal_runs(self):
        """Report two amplicon runs as equal if all their attributes are equal."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        run_1 = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                        configfile = configfile,
                                        active_config = self.active_config)
        run_2 = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                        configfile = configfile,
                                        active_config = self.active_config)
        assert run_1 == run_2

    def test_compare_not_equal_runs(self):
        """Report two amplicon runs as not equal if they differ in one or more attributes."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        run_1 = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                        configfile = configfile,
                                        active_config = self.active_config)
        run_2 = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                        configfile = configfile,
                                        active_config = self.active_config,
                                        sequencing_time = 1.5)
        assert run_1 != run_2

    def test_compare_wrong_type(self):
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        run_1 = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                        configfile = configfile,
                                        active_config = self.active_config)
        OldAmpliconRun = namedtuple("OldAmpliconRun",
                                    ["sequence_dir", "runsheet",
                                     "configfile",
                                     "active_config", "outdir",
                                     "sequencing_time",
                                     "test_run"])
        run_2 = OldAmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                               configfile = configfile,
                               active_config = self.active_config,
                               outdir = pathlib.Path("/path/to/output/"
                                                     "NANO_Amplicon_Y20990101_RUN0001_XYZ-16S"),
                               sequencing_time = 1, test_run = None)
        assert run_1 != run_2
        assert not run_1 == run_2


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
                     "debug": False,
                     "paths": {"output_base_path": "/path/to/output"},
                     "seq_run_duration_hours": 1}
    config_path = pathlib.Path(__file__).parent / "data"/"monitor_run"/"test_config.yaml"
    test_run = monitor_run.AmpliconRun(sequence_dir = pathlib.Path(__file__).parent / "data"
                                                      / "monitor_run" / "miniondir" / "test1",
                                       runsheet = pathlib.Path(__file__).parent / "data"
                                                  / "utilities_test"
                                                  / "test_nanopore_runsheet_16s_only.xlsx",
                                       configfile = config_path, active_config = active_config,
                                       test_run = False)
    test_run.outdir = pathlib.Path(__file__).parent / "data"/ "monitor_run" / "test_outdir"

    def test_timeout(self):
        """Stop monitoring when the timeout has been reached."""
        error_msg = ("No file matching pattern test_summary*.txt found in"
                     f" {self.test_run.sequence_dir / 'rawdata' / 'test_subdir'}")
        with pytest.raises(FileNotFoundError, match=re.escape(error_msg)):
            monitor_run.start_on_file_found(self.test_run, "test_summary*.txt",
                                            watch_interval = 1,
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
        expected_log = ("INFO:launch_run:Found final_summary*.txt in "
                        f"{self.test_run.sequence_dir / 'rawdata' / 'test_subdir'}"
                        " after 0 seconds.")
        with self.assertLogs("launch_run", level = "INFO") as logged:
            test_command = monitor_run.start_on_file_found(self.test_run, "final_summary*.txt",
                                                           dry_run = True, watch_interval = 1,
                                                           watch_timeout = 5).args
            assert expected_log in logged.output
        assert test_command == expected_command

    # we don't want to run the actual nomad command when testing...
    @mock.patch(f'{monitor_run.__name__}.get_pipeline_command',
                wraps = monitor_run.get_pipeline_command)
    def test_success_real_run(self, mock_command):
        """Successfully run the output of the pipeline command builder."""
        mock_command.return_value = ["echo", "Hello"]
        expected_command = ["echo", "Hello"]
        test_command = monitor_run.start_on_file_found(self.test_run, "final_summary*.txt",
                                                       dry_run = False, watch_interval = 1,
                                                       watch_timeout = 5).args
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
                     "debug": False,
                     "paths": {"output_base_path": "/path/to/output"},
                     "seq_run_duration_hours": 1
                     }
    config_path = pathlib.Path(__file__).parent / "data"/"monitor_run"/"test_config.yaml"

    def test_get_routine_command(self):
        indir = pathlib.Path("path/to/indir")
        outdir = pathlib.Path("path/to/outdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        seq_run = monitor_run.AmpliconRun(sequence_dir = indir, runsheet = runsheet,
                                          configfile = self.config_path,
                                          active_config = self.active_config)
        seq_run.outdir = outdir
        expected_command = ["nomad", "job", "dispatch",
                            "-meta", "indir=path/to/indir",
                            "-meta", "outdir=path/to/outdir",
                            "-meta", f"runsheet={runsheet}",
                            "16s-snake-emu-prod", str(self.config_path)]
        test_command = monitor_run.get_pipeline_command(seq_run)
        assert expected_command == test_command

    def test_get_test_command(self):
        """Test that the staging version is run if debug is specified (overriding default)."""
        indir = pathlib.Path("path/to/indir")
        outdir = pathlib.Path("path/to/outdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        expected_command = ["nomad", "job", "dispatch",
                            "-meta", "indir=path/to/indir",
                            "-meta", "outdir=path/to/outdir",
                            "-meta", f"runsheet={runsheet}",
                            "16s-snake-emu-staging", str(self.config_path)]
        seq_run = monitor_run.AmpliconRun(sequence_dir = indir, runsheet = runsheet,
                                          configfile = self.config_path,
                                          active_config = self.active_config, test_run = True)
        seq_run.outdir = outdir
        test_command = monitor_run.get_pipeline_command(seq_run)
        assert expected_command == test_command

    def test_override_test_command(self):
        """Test that a config specifying debug mode can be overridden by the debug flag."""
        indir = pathlib.Path("path/to/indir")
        outdir = pathlib.Path("path/to/outdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        test_configfile = pathlib.Path("path/to/config")
        test_config = {"debug": True, "amplicon_type": "16S",
                       "paths": {"output_base_path": "/path/to/output"},
                       "seq_run_duration_hours": 1}
        expected_command = ["nomad", "job", "dispatch",
                            "-meta", "indir=path/to/indir",
                            "-meta", "outdir=path/to/outdir",
                            "-meta", f"runsheet={runsheet}",
                            "16s-snake-emu-prod", str(test_configfile)]
        seq_run = monitor_run.AmpliconRun(sequence_dir = indir, runsheet = runsheet,
                                          configfile = test_configfile, active_config = test_config,
                                          test_run = False)
        seq_run.outdir = outdir
        test_command = monitor_run.get_pipeline_command(seq_run)
        assert expected_command == test_command
