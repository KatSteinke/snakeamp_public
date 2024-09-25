import copy
import logging
import pathlib
import re
import unittest

from collections import namedtuple
from datetime import timedelta
from unittest import mock

import pytest

import monitor_run


class TestCreateAmpliconRun(unittest.TestCase):
    active_config = {
        "input_names": pathlib.Path(__file__).parent / "data" / "input_names" / "input_da_old_lis.yaml",
        "sample_number_settings": {"sample_number_format":
                                       r'([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})-\d?',
                                                "sample_numbers_in": "letter",
                                                "sample_numbers_out": "letter",
                                                "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                "number_to_letter": {"70": "P",
                                                                     "30": "B",
                                                                     "10": "D",
                                                                     "50": "T"},
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
            monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                    configfile = configfile, active_config = self.active_config,
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

    def test_set_snake_flags(self):
        """Set flags for Snakemake."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        snake_flags = ["-n"]
        amplicon_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config,
                                               snake_flags = snake_flags)
        assert amplicon_run.snake_flags == snake_flags

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
                         f"{timedelta(hours=self.active_config['seq_run_duration_hours'])},"
                         f"snake_flags=None)")
        amplicon_repr = amplicon_run.__repr__()
        assert amplicon_repr == expected_repr

    def test_print_snake_flags(self):
        """Print the run's attributes including Snakemake flags."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        expected_outdir = pathlib.Path("/path/to/output/NANO_Amplicon_Y20990101_RUN0001_XYZ-16S")
        snake_flags = ["-n"]
        amplicon_run = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                               configfile = configfile,
                                               active_config = self.active_config,
                                               snake_flags = snake_flags)
        expected_repr = (f"AmpliconRun(sequence_dir={expected_indir}, runsheet={runsheet}, "
                         f"configfile={configfile}, active_config={self.active_config}, "
                         f"outdir={expected_outdir},"
                         " sequencing_time="
                         f"{timedelta(hours=self.active_config['seq_run_duration_hours'])},"
                         f"snake_flags={snake_flags})")
        amplicon_repr = amplicon_run.__repr__()
        assert amplicon_repr == expected_repr

    def test_compare_equal_runs(self):
        """Report two amplicon runs as equal if all their attributes are equal."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        snake_flags = ["-n", "--rerun-incomplete"]
        run_1 = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                        configfile = configfile, active_config = self.active_config,
                                        snake_flags = snake_flags)
        run_2 = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                        configfile = configfile, active_config = self.active_config,
                                        snake_flags = snake_flags)
        assert run_1 == run_2
        run_3 = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                        configfile = configfile, active_config = self.active_config)
        run_4 = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                        configfile = configfile, active_config = self.active_config)
        assert run_3 == run_4

    def test_compare_not_equal_runs(self):
        """Report two amplicon runs as not equal if they differ in one or more attributes."""
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        snake_flags = ["-n", "--rerun-incomplete"]
        run_1 = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                        configfile = configfile, active_config = self.active_config,
                                        snake_flags = snake_flags)
        run_2 = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                        configfile = configfile, active_config = self.active_config)
        assert run_1 != run_2

    def test_compare_wrong_type(self):
        expected_indir = pathlib.Path("path/to/indir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        configfile = pathlib.Path("path/to/config")
        run_1 = monitor_run.AmpliconRun(sequence_dir = expected_indir, runsheet = runsheet,
                                        configfile = configfile, active_config = self.active_config)
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





class TestGetNomadCommand(unittest.TestCase):
    active_config = {"input_names": pathlib.Path(__file__).parent / "data" / "input_names"
                                    / "input_da_old_lis.yaml",
                     "sample_number_settings": {"sample_number_format":
                                                    r'([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})-\d?',
                                                "sample_numbers_in": "letter",
                                                "sample_numbers_out": "letter",
                                                "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                "number_to_letter": {"70": "P",
                                                                     "30": "B",
                                                                     "10": "D",
                                                                     "50": "T"},

                                                "negative_control": 'NegK[a-zA-Z0-9]*',
                                                "positive_control": {}},
                     "barcode_format": "NB[0-9]{2}",  # format of barcodes in runsheet
                     "barcode_prefix": "NB",
                     "amplicon_type": "16S",
                     "debug": False,
                     "run_on": "nomad",
                     "paths": {"output_base_path": "/path/to/output"},
                     "seq_run_duration_hours": 1
                     }
    config_path = pathlib.Path(__file__).parent / "data"/"monitor_run"/"test_config.yaml"

    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

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
        test_command = monitor_run.get_nomad_command(seq_run)
        assert expected_command == test_command

    def test_warn_unused_flags(self):
        """Warn if run contains snake_flags (unused in the Nomad setup)."""
        indir = pathlib.Path("path/to/indir")
        outdir = pathlib.Path("path/to/outdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        seq_run = monitor_run.AmpliconRun(sequence_dir = indir, runsheet = runsheet,
                                          configfile = self.config_path,
                                          active_config = self.active_config)
        seq_run.outdir = outdir
        seq_run.snake_flags = ["-n"]
        expected_command = ["nomad", "job", "dispatch",
                            "-meta", "indir=path/to/indir",
                            "-meta", "outdir=path/to/outdir",
                            "-meta", f"runsheet={runsheet}",
                            "16s-snake-emu-prod", str(self.config_path)]
        warn_msg = ("Flags to pass to Snakemake are not used when running the pipeline "
                    "in Nomad mode.")
        with self._caplog.at_level(logging.WARNING, logger = "launch_run"):
            test_command = monitor_run.get_nomad_command(seq_run)
            assert ("launch_run", logging.WARNING, warn_msg) in self._caplog.record_tuples
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
        test_command = monitor_run.get_nomad_command(seq_run)
        assert expected_command == test_command

    def test_get_level_from_config(self):
        """Get test level from config file if not specified for the run."""
        indir = pathlib.Path("path/to/indir")
        outdir = pathlib.Path("path/to/outdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        test_configfile = pathlib.Path("path/to/config")
        test_config = {"input_names": pathlib.Path(__file__).parent / "data" / "input_names"
                                      / "input_da_old_lis.yaml",
                       "debug": True, "amplicon_type": "16S",
                       "paths": {"output_base_path": "/path/to/output"},
                       "seq_run_duration_hours": 1}
        expected_command = ["nomad", "job", "dispatch",
                            "-meta", "indir=path/to/indir",
                            "-meta", "outdir=path/to/outdir",
                            "-meta", f"runsheet={runsheet}",
                            "16s-snake-emu-staging", str(test_configfile)]
        seq_run = monitor_run.AmpliconRun(sequence_dir = indir, runsheet = runsheet,
                                          configfile = test_configfile, active_config = test_config,
                                          test_run = None)
        seq_run.outdir = outdir
        test_command = monitor_run.get_nomad_command(seq_run)
        assert expected_command == test_command

    def test_override_test_command(self):
        """Test that a config specifying debug mode can be overridden by the debug flag."""
        indir = pathlib.Path("path/to/indir")
        outdir = pathlib.Path("path/to/outdir")
        runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                    / "test_nanopore_runsheet.xlsx")
        test_configfile = pathlib.Path("path/to/config")
        test_config = {"input_names": pathlib.Path(__file__).parent / "data" / "input_names"
                                      / "input_da_old_lis.yaml",
                       "debug": True, "amplicon_type": "16S",
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
        test_command = monitor_run.get_nomad_command(seq_run)
        assert expected_command == test_command


class TestGetLocalCommand(unittest.TestCase):
    active_config = {"input_names": pathlib.Path(__file__).parent / "data" / "input_names" / "input_da_old_lis.yaml",
                     "sample_number_settings": {"sample_number_format":
                                                    r'([BDPT]|[1357]0)([0-9]{8}|[0-9]{6})-\d?',
                                                "sample_numbers_in": "letter",
                                                "sample_numbers_out": "letter",
                                                "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                "number_to_letter": {"70": "P",
                                                                     "30": "B",
                                                                     "10": "D",
                                                                     "50": "T"},

                                                "negative_control": 'NegK[a-zA-Z0-9]*',
                                                "positive_control": {}},
                     "barcode_format": "NB[0-9]{2}",  # format of barcodes in runsheet
                     "barcode_prefix": "NB",
                     "amplicon_type": "16S",
                     "debug": False,
                     "run_on": "local",
                     "cores": 8,
                     "paths": {"output_base_path": "/path/to/output"},
                     "seq_run_duration_hours": 1,
                     "conda": {"frontend": "",
                               "prefix": "",
                               "use_conda": False}
                     }
    config_path = pathlib.Path(__file__).parent / "data" / "monitor_run" / "test_config.yaml"
    indir = pathlib.Path("path/to/indir")
    outdir = pathlib.Path("path/to/outdir")
    runsheet = (pathlib.Path(__file__).parent / "data" / "utilities_test"
                / "test_nanopore_runsheet.xlsx")
    seq_run = monitor_run.AmpliconRun(sequence_dir = indir, runsheet = runsheet,
                                      configfile = config_path,
                                      active_config = active_config)
    seq_run.outdir = outdir

    def test_get_routine_command(self):
        """Get the command for running the pipeline locally in Snakemake."""
        expected_command = ["snakemake", "-s", "Snakefile",
                            "--cores", "8",
                            "--keep-going",
                            "--config",
                            f"outdir={self.outdir}",
                            f"rundir={self.indir}",
                            f"runsheet={self.runsheet}",
                            f"config_path={str(self.config_path)}",
                            "--configfile", str(self.config_path)]
        test_command = monitor_run.get_local_command(self.seq_run)
        assert expected_command == test_command

    def test_set_snake_flags(self):
        """Get the command for running the pipeline locally in Snakemake."""
        seq_run = monitor_run.AmpliconRun(sequence_dir = self.indir, runsheet = self.runsheet,
                                      configfile = self.config_path,
                                      active_config = self.active_config, outdir = self.outdir)
        seq_run.snake_flags = ["-n"]
        expected_command = ["snakemake", "-s", "Snakefile",
                            "--cores", "8",
                            "--keep-going",
                            "--config",
                            f"outdir={self.outdir}",
                            f"rundir={self.indir}",
                            f"runsheet={self.runsheet}",
                            f"config_path={str(self.config_path)}",
                            "--configfile", str(self.config_path),
                            "-n"]
        test_command = monitor_run.get_local_command(seq_run)
        assert expected_command == test_command

    def test_add_conda_frontend(self):
        """Add conda frontend if given."""
        test_config = copy.deepcopy(self.active_config)
        test_config["conda"]["frontend"] = "mamba"
        test_config["conda"]["use_conda"] = True
        seq_run = monitor_run.AmpliconRun(sequence_dir = self.indir, runsheet = self.runsheet,
                                          configfile = self.config_path,
                                          active_config = test_config, outdir = self.outdir)
        expected_command = ["snakemake", "-s", "Snakefile",
                            "--cores", "8",
                            "--keep-going",
                            "--config",
                            f"outdir={self.outdir}",
                            f"rundir={self.indir}",
                            f"runsheet={self.runsheet}",
                            f"config_path={str(self.config_path)}",
                            "--configfile", str(self.config_path),
                            "--use-conda",
                            "--conda-frontend", "mamba"]
        test_command = monitor_run.get_local_command(seq_run)
        assert expected_command == test_command

    def test_add_conda_prefix(self):
        """Add conda prefix if given."""
        test_config = copy.deepcopy(self.active_config)
        test_config["conda"]["prefix"] = "/data/conda_prefix"
        test_config["conda"]["use_conda"] = True
        seq_run = monitor_run.AmpliconRun(sequence_dir = self.indir, runsheet = self.runsheet,
                                      configfile = self.config_path,
                                      active_config = test_config, outdir = self.outdir)
        expected_command = ["snakemake", "-s", "Snakefile",
                            "--cores", "8",
                            "--keep-going",
                            "--config",
                            f"outdir={self.outdir}",
                            f"rundir={self.indir}",
                            f"runsheet={self.runsheet}",
                            f"config_path={str(self.config_path)}",
                            "--configfile", str(self.config_path),
                            "--use-conda",
                            "--conda-prefix", "/data/conda_prefix"]
        test_command = monitor_run.get_local_command(seq_run)
        assert expected_command == test_command

    def test_add_clusterprofile(self):
        """Add a cluster profile if given."""
        test_config = copy.deepcopy(self.active_config)
        test_config["profile"] = "test/profile"
        seq_run = monitor_run.AmpliconRun(sequence_dir = self.indir, runsheet = self.runsheet,
                                          configfile = self.config_path,
                                          active_config = test_config, outdir = self.outdir)
        expected_command = ["snakemake", "-s", "Snakefile",
                            "--cores", "8",
                            "--keep-going",
                            "--config",
                            f"outdir={self.outdir}",
                            f"rundir={self.indir}",
                            f"runsheet={self.runsheet}",
                            f"config_path={str(self.config_path)}",
                            "--configfile", str(self.config_path),
                            "--profile", "test/profile"]
        test_command = monitor_run.get_local_command(seq_run)
        assert expected_command == test_command


class TestStartGenericRun(unittest.TestCase):
    active_config = {"input_names": pathlib.Path(__file__).parent / "data" / "input_names" / "input_da_old_lis.yaml",
                     "sample_number_settings": {"sample_number_format":
                                                    r'([BDFPT]|[1357]0)([0-9]{8}|[0-9]{6})-\d?',
                                                "sample_numbers_in": "letter",
                                                "sample_numbers_out": "letter",
                                                "format_in_sheet": r'(?P<sample_type>[BDFPT]|[1357]0|11)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                "format_in_lis": r'(?P<sample_type>[BDFPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                "number_to_letter": {"70": "P",
                                                                     "30": "B",
                                                                     "10": "D",
                                                                     "11": "F",
                                                                     "50": "T"},
                                                "negative_control": 'NegK[a-zA-Z0-9]*',
                                                "positive_control": {}},
                     "barcode_format": "NB[0-9]{2}",  # format of barcodes in runsheet
                     "barcode_prefix": "NB",
                     "amplicon_type": "16S",
                     "debug": False,
                     "run_on": "nomad",
                     "cores": 8,
                     "paths": {"output_base_path": "/path/to/output"},
                     "seq_run_duration_hours": 1,
                     "conda": {"frontend": "",
                               "prefix": "",
                               "use_conda": False}}
    config_path = pathlib.Path(__file__).parent / "data" / "monitor_run" / "test_config.yaml"
    test_run = monitor_run.AmpliconRun(sequence_dir = pathlib.Path(__file__).parent / "data"
                                                      / "monitor_run" / "miniondir" / "test1",
                                       runsheet = pathlib.Path(__file__).parent / "data"
                                                  / "utilities_test"
                                                  / "test_nanopore_runsheet_16s_only.xlsx",
                                       configfile = config_path, active_config = active_config,
                                       test_run = False)
    test_run.outdir = pathlib.Path(__file__).parent / "data" / "monitor_run" / "test_outdir"

    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    def test_fail_invalid_mode(self):
        """Fail if an invalid mode is specified."""
        bad_config = self.active_config.copy()
        bad_config["run_on"] = "slurm"
        test_run = copy.copy(self.test_run)
        test_run.active_config = bad_config
        error_msg = "Invalid run mode slurm"
        with pytest.raises(ValueError, match=error_msg):
            monitor_run.start_run_by_mode(test_run, dry_run = True)

    def test_success_dryrun_nomad(self):
        """Print run command for Nomad."""
        expected_command = ['echo', (f'"nomad job dispatch -meta indir={self.test_run.sequence_dir}'
                                     f' -meta outdir={self.test_run.outdir} '
                                     f'-meta runsheet={self.test_run.runsheet} '
                                     '16s-snake-emu-prod '
                                     f'{str(self.test_run.configfile)}"')]
        start_msg = "Running analysis pipeline"
        nomad_msg = "Starting pipeline in nomad mode."
        with self._caplog.at_level(logging.INFO, logger = "launch_run"):
            test_command = monitor_run.start_run_by_mode(self.test_run,
                                                         dry_run = True).args
            assert ("launch_run", logging.INFO, start_msg) in self._caplog.record_tuples
            assert ("launch_run", logging.INFO, nomad_msg) in self._caplog.record_tuples
        assert expected_command == test_command

    def test_success_dryrun_local(self):
        """Print run command for Snakemake."""
        test_run = copy.copy(self.test_run)
        test_config = self.active_config.copy()
        test_config["run_on"] = "local"
        test_run.active_config = test_config
        expected_command = ['echo', (f'"snakemake -s Snakefile '
                                     f'--cores 8 '
                                     f'--keep-going '
                                     f'--config '
                                     f'outdir={self.test_run.outdir} '
                                     f'rundir={self.test_run.sequence_dir} '
                                     f'runsheet={self.test_run.runsheet} '
                                     f'config_path={str(self.config_path)} '
                                     f'--configfile {str(self.config_path)}"')]
        start_msg = "Running analysis pipeline"
        local_msg = "Starting pipeline in local mode."
        with self._caplog.at_level(logging.INFO, logger = "launch_run"):
            test_command = monitor_run.start_run_by_mode(test_run, dry_run = True).args
            assert ("launch_run", logging.INFO, start_msg) in self._caplog.record_tuples
            assert ("launch_run", logging.INFO, local_msg) in self._caplog.record_tuples
        assert expected_command == test_command

    @mock.patch(f'{monitor_run.__name__}.get_nomad_command',
                wraps = monitor_run.get_nomad_command)
    def test_success_real_run(self, mock_command):
        """Successfully run the output of the pipeline command builder."""
        mock_command.return_value = ["echo", "Hello"]
        expected_command = ["echo", "Hello"]
        test_command = monitor_run.start_run_by_mode(self.test_run).args
        assert test_command == expected_command


class TestWaitForFile(unittest.TestCase):
    active_config = {"input_names": pathlib.Path(__file__).parent / "data" / "input_names"
                                    / "input_da_old_lis.yaml",
                     "sample_number_settings": {"sample_number_format":
                                                    r'([BDFPT]|[1357]0)([0-9]{8}|[0-9]{6})-\d?',
                                                "sample_numbers_in": "letter",
                                                "sample_numbers_out": "letter",
                                                "format_in_sheet": r'(?P<sample_type>[BDPT]|[1357]0)(?P<sample_year>\d{2})(?P<sample_number>\d{6})(?P<bact_number>-\d)?',
                                                "format_in_lis": r'(?P<sample_type>[BDPT])(?P<sample_year>\d{2})(?P<sample_number>\d{6})',
                                                "number_to_letter": {"70": "P",
                                                                     "30": "B",
                                                                     "10": "D",
                                                                     "50": "T"},
                                                "negative_control": 'NegK[a-zA-Z0-9]*',
                                                "positive_control": {}},
                     "barcode_format": "NB[0-9]{2}",  # format of barcodes in runsheet
                     "barcode_prefix": "NB",
                     "amplicon_type": "16S",
                     "debug": False,
                     "cores": 8,
                     "run_on": "nomad",
                     "paths": {"output_base_path": "/path/to/output"},
                     "seq_run_duration_hours": 1,
                     "conda": {"frontend": "",
                               "prefix": "",
                               "use_conda": False}}
    config_path = pathlib.Path(__file__).parent / "data"/"monitor_run"/"test_config.yaml"
    test_run = monitor_run.AmpliconRun(sequence_dir = pathlib.Path(__file__).parent / "data"
                                                      / "monitor_run" / "miniondir" / "test1",
                                       runsheet = pathlib.Path(__file__).parent / "data"
                                                  / "utilities_test"
                                                  / "test_nanopore_runsheet_16s_only.xlsx",
                                       configfile = config_path, active_config = active_config,
                                       test_run = False)
    test_run.outdir = pathlib.Path(__file__).parent / "data"/ "monitor_run" / "test_outdir"


    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    def test_timeout(self):
        """Stop monitoring when the timeout has been reached."""
        error_msg = ("No file matching pattern test_summary*.txt found in"
                     f" {self.test_run.sequence_dir / 'no_sample' / 'test_subdir'}")
        time_log = "Waiting for sequencing to finish..."
        log_error = ("No file matching pattern test_summary*.txt "
                     f"found in {self.test_run.sequence_dir / 'no_sample' / 'test_subdir'} "
                     "after 0.0 hours.")
        with (pytest.raises(FileNotFoundError, match=re.escape(error_msg)),
              self._caplog.at_level(level="INFO", logger="launch_run")):
            monitor_run.start_on_file_found(self.test_run, "test_summary*.txt",
                                            watch_interval = 1,
                                            watch_timeout = 3,
                                            log_interval = 1)
        assert ("launch_run", logging.ERROR, log_error) in self._caplog.record_tuples
        timed_logs = [log_record for log_record in self._caplog.record_tuples
                      if log_record == ("launch_run", logging.INFO, time_log)]
        assert len(timed_logs) == 3

    def test_fail_parent_dir_not_found(self):
        """Don't start monitoring if the parent directory does not exist."""
        test_run = monitor_run.AmpliconRun(sequence_dir = pathlib.Path(__file__).parent / "data"
                                                          / "monitor_run" / "miniondir" / "test2",
                                           runsheet = pathlib.Path(__file__).parent / "data"
                                                      / "utilities_test"
                                                      / "test_nanopore_runsheet_16s_only.xlsx",
                                           configfile = self.config_path,
                                           active_config = self.active_config,
                                           outdir = pathlib.Path(__file__).parent / "data"
                                                    / "monitor_run" / "test_outdir",
                                           test_run = False)
        error_msg = f"Parent directory {test_run.sequence_dir} does not exist."
        with pytest.raises(FileNotFoundError, match = re.escape(error_msg)):
            monitor_run.start_on_file_found(test_run, "final_summary*.txt",
                                            watch_interval = 1,
                                            watch_timeout = 5, dry_run = True)

    def test_success_dryrun(self):
        """Successfully print the pipeline start command in dry run mode."""
        expected_command = ['echo', (f'"nomad job dispatch -meta indir={self.test_run.sequence_dir}'
                                     f' -meta outdir={self.test_run.outdir} '
                                     f'-meta runsheet={self.test_run.runsheet} '
                                     '16s-snake-emu-prod '
                                     f'{str(self.test_run.configfile)}"')]
        expected_log = ("Found final_summary*.txt in "
                        f"{self.test_run.sequence_dir / 'no_sample' / 'test_subdir'}"
                        " after 0 seconds.")
        with self._caplog.at_level(logging.INFO, logger = "launch_run"):
            test_command = monitor_run.start_on_file_found(self.test_run, "final_summary*.txt",
                                                           dry_run = True, watch_interval = 1,
                                                           watch_timeout = 5).args
            assert ("launch_run", logging.INFO, expected_log) in self._caplog.record_tuples
        assert test_command == expected_command

    def test_success_dryrun_local(self):
        """Successfully print a local pipeline start command in dry run mode."""
        test_config = copy.deepcopy(self.active_config)
        test_config["run_on"] = "local"
        test_run = copy.copy(self.test_run)
        test_run.active_config = test_config
        expected_command = ['echo', (f'"snakemake -s Snakefile '
                                     f'--cores 8 '
                                     f'--keep-going '
                                     f'--config '
                                     f'outdir={self.test_run.outdir} '
                                     f'rundir={self.test_run.sequence_dir} '
                                     f'runsheet={self.test_run.runsheet} '
                                     f'config_path={str(self.config_path)} '
                                     f'--configfile {str(self.config_path)}"')]
        expected_log = ("Found final_summary*.txt in "
                        f"{self.test_run.sequence_dir / 'no_sample' / 'test_subdir'}"
                        " after 0 seconds.")
        with self._caplog.at_level(logging.INFO, logger = "launch_run"):
            test_command = monitor_run.start_on_file_found(test_run, "final_summary*.txt",
                                                           dry_run = True, watch_interval = 1,
                                                           watch_timeout = 5).args
            assert ("launch_run", logging.INFO, expected_log) in self._caplog.record_tuples
        assert test_command == expected_command

    # we don't want to run the actual nomad command when testing...
    @mock.patch(f'{monitor_run.__name__}.get_nomad_command',
                wraps = monitor_run.get_nomad_command)
    def test_success_real_run(self, mock_command):
        """Successfully run the output of the pipeline command builder."""
        mock_command.return_value = ["echo", "Hello"]
        expected_command = ["echo", "Hello"]
        test_command = monitor_run.start_on_file_found(self.test_run, "final_summary*.txt",
                                                       dry_run = False, watch_interval = 1,
                                                       watch_timeout = 5).args
        assert test_command == expected_command
