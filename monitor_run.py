"""Monitor the output directory of a specified sequencing run and start 16S analysis pipeline."""

__author__ = "Kat Steinke"

import logging
import pathlib
import subprocess
import time

from typing import Any, Dict, Optional, List

import helpers
import pipeline_config
import version

__version__ = version.__version__

# import parameters
default_config_file = pipeline_config.default_config_file
workflow_config = pipeline_config.WORKFLOW_DEFAULT_CONF

# start logging
logger = logging.getLogger("launch_run")
logger.setLevel(logging.INFO)
console_log = logging.StreamHandler()
console_log.setLevel(logging.INFO)
logger.addHandler(console_log)


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
    """
    def __init__(self, sequence_dir: pathlib.Path, runsheet: pathlib.Path, configfile: pathlib.Path,
                 active_config: Dict[str, Any], outdir: Optional[pathlib.Path] = None,
                 sequencing_time: Optional[float] = None,
                 test_run: Optional[bool] = None):
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
        """
        self.sequence_dir = sequence_dir
        self.runsheet = runsheet
        self.configfile = configfile
        self.active_config = active_config
        self.test_run = test_run
        if outdir:
            self.outdir = outdir
        else:
            self.outdir = (pathlib.Path(active_config['paths']['output_base_path'])
                           / f"{helpers.extract_nanopore_run_name(runsheet)}"
                             f"-{active_config['amplicon_type']}")
        if sequencing_time:
            self.sequencing_time = float(sequencing_time)
        else:
            self.sequencing_time = float(active_config['seq_run_duration_hours'])


def get_pipeline_command(sequencing_run: AmpliconRun) -> List[str]:
    """Generate the command for starting the pipeline.

    Arguments:
        sequencing_run:          the directory containing input files for the pipeline

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
    return nomad_command


# convenience function to allow user to specify duration

def start_on_file_found(run_to_watch: AmpliconRun, pattern_to_watch: str, dry_run: bool = False,
                        watch_interval: int = 300,
                        watch_timeout: int = 600) -> subprocess.CompletedProcess:
    """Start the analysis pipeline when a given file is found in the specified run's sequencing dir.

    Arguments:
        run_to_watch:       the run for which the analysis pipeline should be started
        pattern_to_watch:   the file pattern to watch for
        dry_run:            whether or not to only print the pipeline command
        watch_interval:     the interval (in seconds) in which the script should check for
                            the presence of the file
        watch_timeout:      the timespan (in seconds) to wait for the file


    Returns:
        The process that launches the pipeline
    Raises:
        FileNotFoundError:  if the file is not found before the timeout
    """
    file_found = 0
    time_watching = 0
    run_basedir = run_to_watch.sequence_dir / "rawdata"
    # watch for presence of file
    while not (file_found or time_watching >= watch_timeout):
        file_found = len(list(run_basedir.glob(pattern_to_watch)))
        if not file_found:  # TODO: we can definitely make this flow more nicely
            time_watching += watch_interval
            time.sleep(watch_interval)
    if not file_found:
        raise FileNotFoundError(f"No file matching pattern {pattern_to_watch} "
                                f"found in {run_basedir}.")
    logger.info(f"Found {pattern_to_watch} in {run_basedir} after {time_watching} seconds.")
    # construct nomad command
    nomad_command = get_pipeline_command(run_to_watch)
    # return only string if in test mode
    if dry_run:
        nomad_command_text = " ".join(nomad_command)
        return subprocess.run(["echo", f'"{nomad_command_text}"'])
    return subprocess.run(nomad_command)

# parser to take input? but we need to pass all parameters down
