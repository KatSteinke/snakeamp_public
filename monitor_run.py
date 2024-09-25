"""Monitor the output directory of a specified sequencing run and start 16S analysis pipeline."""

__author__ = "Kat Steinke"

#  Copyright (c) 2026 Kat Steinke
#     This program is distributed under version 3 of the GNU General Public License.
#      You should have received a copy of the GNU General Public License
#        along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import logging
import pathlib
import subprocess
import time

from datetime import timedelta
from typing import Any, Dict, Optional, List

import helpers
import input_names
import pipeline_config
import version

__version__ = version.__version__

# import parameters
default_config_file = pipeline_config.default_config_file
workflow_config = pipeline_config.WORKFLOW_DEFAULT_CONF

# start logging
logger = logging.getLogger("launch_run")


class AmpliconRun:
    """Parameters for an amplicon sequencing run.

    Attributes:
        sequence_dir:       the directory containing input files for the pipeline
        outdir:             the directory to which results should be output
        runsheet:           the runsheet used for the run
        configfile:         the file containing the configuration for the pipeline
        active_config:      the configuration to use for the pipeline
        outdir:             the output directory to which results should be output
                            (inferred from experiment name if not given)
        sequencing_time:    the expected duration of the run in hours
                            (taken from config if not given)
        test_run:           whether to run the pipeline in test mode (overrides config setting)
        snake_flags         optional flags to pass on to Snakemake
    """
    def __init__(self, sequence_dir: pathlib.Path, runsheet: pathlib.Path, configfile: pathlib.Path,
                 active_config: Dict[str, Any], outdir: Optional[pathlib.Path] = None,
                 sequencing_time: Optional[float] = None, test_run: Optional[bool] = None,
                 snake_flags: Optional[List[str]] = None):
        """Initialize an AmpliconRun with the supplied parameters, setting the output dir to one
        based on the experiment name.

        Arguments:
            sequence_dir:       the directory containing input files for the pipeline
            runsheet:           the runsheet used for the run
            configfile:         the file containing the configuration for the pipeline
            active_config:      the configuration to use for the pipeline
            outdir:             the output directory to which results should be output
                                (inferred from experiment name if not given)
            sequencing_time:    the expected duration of the run in hours
                                (taken from config if not given)
            test_run:           whether to run the pipeline in test mode (overrides config setting)
            snake_flags:        flags to pass on to Snakemake
        """
        self.sequence_dir = sequence_dir
        self.runsheet = runsheet
        self.configfile = configfile
        self.active_config = active_config
        self.snake_flags = snake_flags
        if test_run is None:
            self.test_run = active_config["debug"]
        else:
            self.test_run = test_run
        if outdir:
            self.outdir = outdir
        else:
            _runsheet_names, _lis_names = input_names.load_input_from_config(active_config)
            self.outdir = (pathlib.Path(active_config['paths']['output_base_path'])
                           / f"{helpers.extract_nanopore_run_name(runsheet, _runsheet_names)}"
                             f"-{active_config['amplicon_type']}")
        seq_time = float(active_config['seq_run_duration_hours'])
        if sequencing_time:
            seq_time = float(sequencing_time)
        if seq_time < 0:
            raise ValueError("Expected sequencing time must be greater than 0 hours.")
        self.sequencing_time = timedelta(hours=seq_time)

    def __repr__(self):
        return (f"AmpliconRun(sequence_dir={self.sequence_dir}, runsheet={self.runsheet}, "
                f"configfile={self.configfile}, active_config={self.active_config}, "
                f"outdir={self.outdir}, sequencing_time={self.sequencing_time},"
                f"snake_flags={self.snake_flags})")

    def __eq__(self, other):
        if isinstance(other, AmpliconRun):
            return ((self.sequence_dir == other.sequence_dir
                    and self.runsheet == other.runsheet
                    and self.configfile == other.configfile
                    and self.active_config == other.active_config
                    and self.outdir == other.outdir
                    and self.sequencing_time == other.sequencing_time
                    and self.test_run == other.test_run
                     and self.snake_flags == other.snake_flags))
        return False

    def __ne__(self, other):
        if isinstance(other, AmpliconRun):
            return (not(self.sequence_dir == other.sequence_dir
                    and self.runsheet == other.runsheet
                    and self.configfile == other.configfile
                    and self.active_config == other.active_config
                    and self.outdir == other.outdir
                    and self.sequencing_time == other.sequencing_time
                    and self.test_run == other.test_run
                        and self.snake_flags == other.snake_flags))
        return True


def get_nomad_command(sequencing_run: AmpliconRun) -> List[str]:
    """Generate the command for starting the pipeline on Nomad.

    Arguments:
        sequencing_run: the sequencing run for which to run the pipeline

    Returns:
        The nomad command to start the pipeline
    """
    if sequencing_run.test_run is None:
        run_as_debug = sequencing_run.active_config["debug"]
    else:
        run_as_debug = sequencing_run.test_run
    nomad_job = "16s-snake-emu-staging" if run_as_debug else "16s-snake-emu-prod"
    nomad_command = ["nomad", "job", "dispatch",
                     "-meta", f"indir={sequencing_run.sequence_dir}",
                     "-meta", f"outdir={sequencing_run.outdir}",
                     "-meta", f"runsheet={sequencing_run.runsheet}",
                     nomad_job,
                     str(sequencing_run.configfile)]

    # might be interesting to pass more flags to Snakemake on Nomad but not implemented for now
    if sequencing_run.snake_flags:
        logger.warning("Flags to pass to Snakemake are not used"
                       " when running the pipeline in Nomad mode.")

    return nomad_command


def get_local_command(sequencing_run: AmpliconRun) -> List[str]:
    """Generate the command for starting the pipeline locally.

    Arguments:
        sequencing_run: the sequencing run for which to run the pipeline

    Returns:
        The snakemake command to start the pipeline.
    """
    snakemake_command = ["snakemake", "-s", "Snakefile",
                        "--cores", str(sequencing_run.active_config["cores"]),
                        "--keep-going",
                        "--config",
                        f"outdir={sequencing_run.outdir}",
                        f"rundir={sequencing_run.sequence_dir}",
                        f"runsheet={sequencing_run.runsheet}",
                        f"config_path={str(sequencing_run.configfile)}",
                        "--configfile", str(sequencing_run.configfile)]
    profile_settings = sequencing_run.active_config.get("profile")
    if profile_settings:
        snakemake_command.extend(["--profile", profile_settings])
    conda_extras = sequencing_run.active_config["conda"]
    if conda_extras["use_conda"]:
        snakemake_command.append("--use-conda")
    if conda_extras["frontend"]:
        snakemake_command.extend(["--conda-frontend", conda_extras["frontend"]])
    if conda_extras["prefix"]:
        snakemake_command.extend(["--conda-prefix", conda_extras["prefix"]])
    if sequencing_run.snake_flags:
        snakemake_command.extend(sequencing_run.snake_flags)
    return snakemake_command


def start_run_by_mode(sequencing_run: AmpliconRun,
                      dry_run: bool = False) -> subprocess.CompletedProcess:
    """Start the analysis pipeline for a given amplicon sequencing run,
     on Nomad or locally as specified in the run's config.

    Arguments:
        sequencing_run: the run for which the analysis pipeline should be started
        dry_run:        whether or not to only print the pipeline command

    Returns:
        The process that launches the pipeline

    Raises:
        ValueError: if the run mode is invalid (not local or nomad)
    """
    pipeline_mode = sequencing_run.active_config["run_on"]
    if pipeline_mode == "nomad":
        pipeline_command = get_nomad_command(sequencing_run)
    elif pipeline_mode == "local":
        pipeline_command = get_local_command(sequencing_run)
    else:
        raise ValueError(f"Invalid run mode {pipeline_mode}")
    logger.info("Running analysis pipeline")
    logger.info(f"Starting pipeline in {pipeline_mode} mode.")
    # return only string if in test mode
    if dry_run:
        pipeline_command_text = " ".join(pipeline_command)
        return subprocess.run(["echo", f'"{pipeline_command_text}"'])
    return subprocess.run(pipeline_command)


# allow user to specify duration
def start_on_file_found(run_to_watch: AmpliconRun, pattern_to_watch: str, dry_run: bool = False,
                        watch_interval: int = 300, watch_timeout: int = 600,
                        log_interval: int = 300) -> subprocess.CompletedProcess:
    """Start the analysis pipeline when a given file is found in the specified run's sequencing dir.

    Arguments:
        run_to_watch:       the run for which the analysis pipeline should be started
        pattern_to_watch:   the file pattern to watch for
        dry_run:            whether or not to only print the pipeline command
        watch_interval:     the interval (in seconds) in which the script should check for
                            the presence of the file
        watch_timeout:      the timespan (in seconds) to wait for the file
        log_interval:       the interval (in seconds) in which the script should log activity



    Returns:
        The process that launches the pipeline
    Raises:
        FileNotFoundError:  if the file is not found before the timeout
    """
    seconds_per_hour = 3600

    file_found = 0
    time_watching = 0
    time_since_log = 0
    if not run_to_watch.sequence_dir.exists():
        raise FileNotFoundError(f"Parent directory {run_to_watch.sequence_dir} does not exist.")
    fastq_parent_dir = helpers.get_fastq_pass_parent(run_to_watch.sequence_dir)
    # watch for presence of file
    while not (file_found or time_watching >= watch_timeout):
        # keep a log so we know it's alive
        time_since_log += watch_interval
        if time_since_log >= log_interval:
            logger.info("Waiting for sequencing to finish...")
            time_since_log = 0
        file_found = len(list(fastq_parent_dir.glob(pattern_to_watch)))
        if not file_found:  # TODO: we can definitely make this flow more nicely
            time_watching += watch_interval
            time.sleep(watch_interval)
    if not file_found:
        timeout_hours = watch_timeout / seconds_per_hour
        logger.error(f"No file matching pattern {pattern_to_watch} found in {fastq_parent_dir} "
                     f"after {timeout_hours:.1f} hours.")
        raise FileNotFoundError(f"No file matching pattern {pattern_to_watch} "
                                f"found in {fastq_parent_dir}.")
    logger.info(f"Found {pattern_to_watch} in {fastq_parent_dir} after {time_watching} seconds.")
    # get command
    return start_run_by_mode(run_to_watch, dry_run = dry_run)
