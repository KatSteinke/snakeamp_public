"""Wrapper script to start the pipeline from the terminal."""

__author__ = "Kat Steinke"

import argparse
import logging
import os
import pathlib
import re
import readline
import subprocess
import sys

from argparse import ArgumentParser
from datetime import timedelta
from typing import Any, Dict, List, Optional

import pandas as pd
import yaml

import check_runsheet
import helpers
import monitor_run
import pipeline_config
import set_log
import version

__version__ = version.__version__


# import parameters
DEFAULT_CONFIG_FILE = pipeline_config.default_config_file
WORKFLOW_CONFIG = pipeline_config.WORKFLOW_DEFAULT_CONF

# start logging
# TODO: capture warnings etc
logger = logging.getLogger("amplicon_nanopore")


# get base dir - see SARS script
# handle bad filenames separately
class BadPathError(Exception):
    """Exception raised when a path to be created contains illegal characters or reserved
    filenames. Workaround to be able to distinguish from incorrect user input in classic mode.
    """
    pass


# TODO: how long does it take for the run dir to be created?

# read runsheet
def process_runsheet(runsheet_path: pathlib.Path,
                     active_config: Dict[str, Any] = WORKFLOW_CONFIG) -> pd.DataFrame:
    """Read a runsheet, filter it down to the samples for which the analysis specified in config
    should be performed and check the format.

    Arguments:
        runsheet_path:  the path to the runsheet to process
        active_config:  the config in use

    Returns:
        The filtered, checked runsheet.

    Raises:
        KeyError:   if there are no samples for which the analysis should be performed
    """
    run_data = pd.read_excel(runsheet_path, usecols = "A:D", skiprows = 3,  # don't check CP for now
                             dtype = {"Prøvenummer": str, "Eluat nr.": str})
    run_data = run_data.dropna(subset = "Prøvenummer")
    # filter down to correct analysis
    run_data = run_data[run_data["Analyse"] == active_config["amplicon_type"]]
    if run_data.empty:
        raise KeyError(f"No samples with amplicon type {active_config['amplicon_type']} "
                       "found in runsheet")
    check_runsheet.check_sheet_format(run_data, check_barcodes = True,
                                      active_config = active_config)
    return run_data

# check runsheet against LIS
# complain if a sample isn't in the LIS report and not a recorded control


# get experiment name and infer output dir

# ensure our output dir is clean
def get_existing_path(path_to_check: pathlib.Path) -> pathlib.Path:
    """Recursively check if path exists, else go down one level until an existing path is found.
    Adapted from https://stackoverflow.com/a/39489505
    Arguments:
        path_to_check:  Path whose components should be checked
    Returns:
        The existing parts of the path
    """
    if path_to_check.exists():
        return path_to_check
    return get_existing_path(path_to_check.parent)


def get_clean_outdir(outdir_path: pathlib.Path) -> pathlib.Path:
    """Check which parts of a path already exist and sanitize the new ones by removing spaces
    and special characters if needed.
    Arguments:
        outdir_path:    Path to check for spaces and special characters
    Returns:
         The sanitized version of the path
    """
    illegal_in_windows = r'[<>:"|?*]'
    # "magic" filenames in Windows, should not be used
    #device_names = {"CON", "PRN", "AUX", "NUL", "COM0", "COM1", "COM2", "COM3", "COM4", "COM5",
    #                "COM6", "COM7", "COM8", "COM9", "LPT0", "LPT1", "LPT2", "LPT3", "LPT4",
    #                "LPT5", "LPT6", "LPT7", "LPT8", "LPT9"}
    # find out until which point path exists
    existing_path = get_existing_path(outdir_path)
    # for anything up to that, if it contains a "bad" character or a space, fail immediately
    if " " in str(existing_path):
        raise BadPathError("The path you are trying to save results to contains a space "
                           "in an existing folder's name. This can break the pipeline. "
                           "\nAborting....")
    if re.search(illegal_in_windows, str(existing_path)):
        raise BadPathError("The path you are trying to save results to contains a character that "
                           "can't be used in Windows in an existing folder's name. "
                           "This can break the pipeline. "
                           "\nAborting....")
    if existing_path == outdir_path:
        return existing_path
    # get rest of the path: relative to existing path
    new_path = outdir_path.relative_to(existing_path)
    # for the new part of the filename:
    # otherwise remove the "illegal" characters and substitute spaces with underscores
    # TODO: any way to use pathlib for this?
    plain_path = str(new_path)
    plain_path = plain_path.replace(" ", "_")
    plain_path = re.sub(illegal_in_windows, "", plain_path)
    cleaned_path = existing_path / plain_path
    # if the path contains a device name, stop and complain
    # if device_names.intersection(set(new_path.parts)):
    if pathlib.PureWindowsPath(cleaned_path).is_reserved():
        raise BadPathError("The path you are trying to save results to contains a name that is "
                           "reserved in Windows. Cannot create this path. \n"
                           "Aborting....")
    return cleaned_path


# get sequencing time from manual input - TODO: just give it as hours and let the class handle the delta
def ask_seq_time(default_seq_time: float) -> float:
    """Ask the user whether sequencing time is correct and get changed sequencing time if needed

    Arguments:
        default_seq_time:   the default sequencing time specified in the config

    Returns:
        The sequencing timespan for the run
    """
    seq_time_accept = input("Expecting sequencing to be finished after"
                            f" {default_seq_time} hours. "
                            "Is this correct? [y/n]")
    # if they just accept we're done
    if seq_time_accept == "y":
        return default_seq_time
    if seq_time_accept == "n":
        sequencing_time = input("Type how many hours the sequencing run is expected to last"
                                " (e.g. 2 if you set it to 2 hours) and press enter: ")
        try:
            sequencing_time = float(sequencing_time)
        except ValueError as value_err:
            raise ValueError(f"{sequencing_time} is not a valid sequencing time. "
                             "Sequencing time must be entered as numbers "
                             "(e.g. 8 for eight hours or 0.5 for half an hour).") from value_err
        if sequencing_time < 0:
            raise ValueError("Expected sequencing time must be greater than 0 hours.")
        return sequencing_time
    # we should not reach this with valid input
    raise ValueError("Sequencing time not entered. Aborting")


# get output dir from manual input
def ask_output_dir(output_dir_path: pathlib.Path) -> pathlib.Path:
    """Ask the user to confirm the default output directory or define a new one.

    Arguments:
        output_dir_path:    the suggested output directory

    Returns:
        The output directory to use for the run

    Raises:
        BadPathError:   if the path that was entered would break on a Windows file system
        ValueError:     if an invalid response is entered
    """
    try:
        # check if this would break anything in Windows or contains spaces
        output_dir_path = get_clean_outdir(output_dir_path)
        # if it's good, ask user for confirmation
        target_accept = input(f"Do you accept {str(output_dir_path)} as target folder [y/n]")
        if target_accept == "n":
            output_dir_path = pathlib.Path(
                input("Type full path or name of target folder and press "
                      "enter: ").strip().strip("'"))
        elif target_accept == "y":  # TODO: can we handle this more nicely?
            pass
        else:
            raise ValueError("Output folder not entered. Aborting")
    except BadPathError:
        logger.warning("The default target folder contains characters that can break the pipeline.")
        output_dir_path = pathlib.Path(input("Type full path to new target folder "
                                             "and press enter: ").strip().strip("'"))
    # if something still is broken now, we yell at the user and fail
    output_dir_path = get_clean_outdir(output_dir_path)
    logger.info(f"Saving results to {output_dir_path}")
    return output_dir_path


def initialize_classic_run(active_config: Dict[str, Any],
                           configfile: pathlib.Path) -> monitor_run.AmpliconRun:
    """Initialize an amplicon sequencing run in classic mode (asking for user input).

    Arguments:
        active_config:  the configuration to use
        configfile:     the configuration file for the configuration given

    Returns:
        An AmpliconRun with the input directory and runsheet specified by the user,
        with default output directory for now. Sequencing time can be specified, otherwise defaults
        to the time set in config.

    """
    # greet the user - TODO: simplified logging here or generally?
    logger.info(f"### Nanopore {active_config['amplicon_type']} analysis\n"
                "# Setup analysis -------------------------------")
    run_dir = pathlib.Path(input("Type full path or name of Nanopore "
                                 "sequencing folder and press enter: ").strip().strip("'"))
    run_dir = helpers.get_fastq_pass_parent(run_dir)
    sequencing_time = ask_seq_time(active_config["seq_run_duration_hours"])
    run_sheet = pathlib.Path(input("Output directory will be based on experiment name."
                                   "\n"
                                   "Enter path to runsheet: ").strip().strip("'")).resolve()
    basic_run = monitor_run.AmpliconRun(sequence_dir = run_dir, runsheet = run_sheet,
                                        configfile = configfile,
                                        active_config = active_config,
                                        sequencing_time = sequencing_time)
    outdir = ask_output_dir(basic_run.outdir)
    basic_run.outdir = outdir
    return basic_run


def initialize_commandline_run(start_args: argparse.Namespace,
                               active_config: Dict[str, Any] = WORKFLOW_CONFIG,
                               config_file: pathlib.Path = DEFAULT_CONFIG_FILE) \
        -> monitor_run.AmpliconRun:
    """Initialize a run from commandline arguments.

    Arguments:
        start_args:     the arguments the pipeline was started with
        active_config:  the configuration to use
        config_file:    the path to the config file used

    Returns:
        An AmpliconRun with all relevant features taken from commandline arguments.
    """
    if not start_args.runsheet:
        raise ValueError("Runsheet not specified.")
    if not start_args.rundir:
        raise ValueError("Nanopore run directory not specified.")
    # fetch LIS stuff if needed
    seqtime = active_config["seq_run_duration_hours"]
    test_run = active_config["debug"]
    run_sheet = pathlib.Path(start_args.runsheet).resolve()
    run_dir = helpers.get_fastq_pass_parent(pathlib.Path(start_args.rundir))
    if start_args.run_time:
        seqtime = start_args.run_time
    if start_args.test_run:
        test_run = True
    current_amplicon_run = monitor_run.AmpliconRun(sequence_dir = run_dir, runsheet = run_sheet,
                                                   sequencing_time = seqtime,
                                                   configfile = config_file,
                                                   active_config = active_config,
                                                   test_run = test_run)
    if start_args.outdir:  # can only be given in commandline mode
        current_amplicon_run.outdir = get_clean_outdir(pathlib.Path(start_args.outdir))
    return current_amplicon_run


def set_up_output(outdir: pathlib.Path, continue_run: bool = False) -> None:
    """Set up output and log directories if the output directory does not exist already.

    Arguments:
        outdir:         the desired output directory
        continue_run:   whether this is a continued run from an existing output directory

    Raises:
        FileExistsError:    if the output directory already exists and this should be a new run

    """
    if not outdir.exists():
        outdir.mkdir(parents = True)
    else:
        if not continue_run:
            raise FileExistsError("The desired output directory already exists.")
        # if we're continuing, confirm this to the user
        logger.info("Output directory already exists. Continuing run...")
    # create dir for snakemake logs
    if not (outdir / "logs").exists():
        (outdir / "logs").mkdir()

# TODO move this out
def check_if_classic_mode(start_args: argparse.Namespace, parent_parser: ArgumentParser,
                          classic_mode_args: Optional[List[str]] = None) -> bool:
    """Check if the pipeline was called in classic mode or commandline mode. In classic mode,
    all arguments are assumed to be default with the exception of those given in classic_mode_args.

    Arguments:
        start_args:         the arguments the pipeline was called with
        parent_parser:      the parser used to parse the args
        classic_mode_args:  names of the args allowed in classic mode, if any

    Returns:
        True if the pipeline was called with all default args
        (except args which may be set in classic mode);
        False otherwise.

    Raises:
        KeyError:   if the start args or args permitted in classic mode contain args not found in
                    the parent parser
    """
    args_with_vals = vars(start_args)
    parser_actions = parent_parser._actions
    parser_args = {arg.dest for arg in parser_actions}
    excess_args = set(args_with_vals.keys()) - parser_args
    if excess_args:
        error_msg = helpers.PrettyKeyErrorMessage(f"Argument(s) {sorted(list(excess_args))}"
                                                  " are not valid arguments for the "
                                                  "specified parser.")
        raise KeyError(error_msg)
    # we may have args that can be set in classic mode - if these are correct...
    if classic_mode_args:
        extra_classic_args = set(classic_mode_args) - parser_args
        if extra_classic_args:
            error_msg = helpers.PrettyKeyErrorMessage("Argument(s) "
                                                      f"{sorted(list(extra_classic_args))}"
                                                      " are not valid arguments for the "
                                                      "specified parser.")
            raise KeyError(error_msg)
        # ...remove them from what we check
        args_with_vals = {key: value for key, value in args_with_vals.items()
                          if key not in classic_mode_args}

    # start with the easiest case: everything is None/False
    # - we can only have this if no default is True
    if not any(arg.default is True for arg in parser_actions):
        if not any(args_with_vals.values()):
            return True
    # however, even if there are values, these might be defaults
    # https://stackoverflow.com/a/44543594/15704972 for getting parser defaults
    parent_defaults = {key: parent_parser.get_default(key) for key in args_with_vals}
    all_default = args_with_vals == parent_defaults
    return all_default


def run_pipeline(start_args: List[str]) -> subprocess.CompletedProcess:
    """Run the pipeline.

    Arguments:
        start_args: The arguments to start the pipeline with

    Returns:
        The subprocess that runs the pipeline.
    """
    arg_parser = ArgumentParser(description = "Run the Nanopore amplicon analysis pipeline")
    arg_parser.add_argument("--rundir", help="Full path or name of sequencing folder")
    arg_parser.add_argument("--runsheet", help="Path to runsheet")
    arg_parser.add_argument("--outdir",
                            help="Path to output directory "
                                 "(default: "
                                 f"{WORKFLOW_CONFIG['paths']['output_base_path']}/[name of rundir])")
    arg_parser.add_argument("--workflow_config_file",
                            help="Config file for run (overrides default config given in script, "
                                 "can be overridden by commandline options)")
    arg_parser.add_argument("--continue_pipeline", action="store_true",
                            help="Continue pipeline after interruption")
    arg_parser.add_argument("--test_run", action="store_true",
                            help="Start a run with the development version of the pipeline")
    arg_parser.add_argument("--dry_run", action="store_true",
                            help="Determine the pipeline start command,"
                                 " set up required dirs and exit")
    arg_parser.add_argument("--run_time", action="store", type=float,
                            help = "Expected sequencing time in hours "
                                   "(will wait for the sequencing run for another hour after this;"
                                   f" default: {WORKFLOW_CONFIG['seq_run_duration_hours']})")
    arg_parser.add_argument("--logfile", action = "store",
                            help = "Logfile to store analysis start commands (default: "
                                   f"{WORKFLOW_CONFIG['paths']['output_base_path']}/"
                                   "[experiment_name]/logs/start_pipeline.log",
                            default = None)
    args = arg_parser.parse_args(start_args)
    # initialize root logger
    pipeline_logger = set_log.get_stream_log()
    # we'll save to file later, but some things will be logged before we know where to save them to
    # -> store them in the meantime
    log_store = set_log.RecordsListHandler()
    pipeline_logger.addHandler(log_store)
    # initialize defaults
    active_config = WORKFLOW_CONFIG
    config_file = DEFAULT_CONFIG_FILE
    # classic mode is allowed to have a config file
    # load config if present - we need to do this early
    if args.workflow_config_file:
        config_file = pathlib.Path(args.workflow_config_file).resolve()
        with open(config_file, "r", encoding = "utf-8") as read_config:
            active_config = yaml.safe_load(read_config)
    manual_mode = check_if_classic_mode(args, arg_parser,
                                        classic_mode_args = ["workflow_config_file"])


    # load debug settings from config (either the one we loaded or the default)
    debug_run = WORKFLOW_CONFIG["debug"]
    # load sequencing time settings
    seq_time = WORKFLOW_CONFIG["seq_run_duration_hours"]
    seq_run_fudge_factor = timedelta(hours = 1)
    if manual_mode:
        # lots of typing, so allow tab completion of paths
        readline.set_completer_delims('\t\n=')
        readline.parse_and_bind("tab: complete")
        # set logger to something quiet
        plain_messages = logging.Formatter("%(message)s")
        pipeline_logger.handlers[0].setFormatter(plain_messages)
        current_run = initialize_classic_run(active_config, config_file)
    else:
        current_run = initialize_commandline_run(args, active_config, config_file)

    # check runsheet
    runsheet_data = process_runsheet(current_run.runsheet, WORKFLOW_CONFIG)
    # set up use of LIS features if enabled - TODO: do we only use them for the runsheet check?
    if active_config["lab_info_system"]["use_lis_features"]:
        lis_report = active_config["lab_info_system"]["lis_report"]
        check_runsheet.check_against_lis(runsheet_data, lis_report, active_config = WORKFLOW_CONFIG)
    else:
        lis_report = None

    # manual mode has set up the output dir, commandline may still have to
    if args.outdir:  # can only be given in commandline mode
        current_run.outdir = get_clean_outdir(pathlib.Path(args.outdir))

    # we can check for whether this is a test and/or a continued pipeline here - in manual mode the
    # arguments will be false, so it defaults to pipeline defaults
    if args.test_run:
        debug_run = True
    # now it's set for sure we can set it on the AmpliconRun
    current_run.test_run = debug_run
    continue_pipeline = args.continue_pipeline
    # we'll have to handle creating our folders ourselves - catch duplicate dirs here!
    # TODO can we move this up at all?
    #  Right now if we do the pipeline will balk if you've fixed a broken runsheet
    #  and try to run with default params in classic mode again
    set_up_output(outdir = current_run.outdir, continue_run = continue_pipeline)
    if args.logfile:
        log_file = pathlib.Path(args.logfile)
    else:
        log_file = current_run.outdir / "logs" / "start_pipeline.log"
    file_handle = logging.FileHandler(log_file)
    # file logs ought to be a bit more detailed
    file_format = logging.Formatter(fmt = '%(asctime)s - %(levelname)s: %(module)s: %(message)s')
    file_handle.setFormatter(file_format)
    # properly attach the handler
    pipeline_logger.addHandler(file_handle)
    # add all of our old messages to the file
    set_log.handle_bulk_logs(file_handle, log_store.record_log)
    # ... and remove our temporary one
    pipeline_logger.removeHandler(log_store)
    log_store.close()
    # start pipeline (in Docker container)
    logger.info("Pipeline is now waiting for sequencing to finish...")
    # set up sequencing run
    # calculate waiting time and check interval
    total_time = current_run.sequencing_time + seq_run_fudge_factor
    check_interval = active_config["check_interval_seconds"]
    # wait and start - TODO: give pattern more nicely?
    # detach here - keep start log
    if os.fork():
        sys.exit()
    analysis_run = monitor_run.start_on_file_found(current_run, "final_summary*.txt",
                                                   dry_run = args.dry_run,
                                                   watch_timeout = int(total_time.total_seconds()),
                                                   watch_interval = check_interval)  # TODO: add logging interval
    logger.info(f"Started pipeline with command {' '.join(analysis_run.args)}")
    # clean up the remaining handlers
    set_log.clean_up_handlers(pipeline_logger)
    return analysis_run


# get the command to run the pipeline
if __name__ == "__main__":
    run_pipeline(sys.argv[1:])
