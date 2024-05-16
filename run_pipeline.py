"""Wrapper script to start the pipeline from the terminal."""

__author__ = "Kat Steinke"

import logging
import os
import pathlib
import re
import readline
import sys

from argparse import ArgumentParser
from datetime import timedelta
from typing import Any, Dict

import pandas as pd
import yaml

import check_runsheet
import helpers
import monitor_run
import pipeline_config
import version

__version__ = version.__version__


# import parameters
default_config_file = pipeline_config.default_config_file
workflow_config = pipeline_config.WORKFLOW_DEFAULT_CONF

# start logging
logger = logging.getLogger("amplicon_nanopore")
logger.setLevel(logging.INFO)
console_log = logging.StreamHandler()
console_log.setLevel(logging.INFO)
logger.addHandler(console_log)


# get base dir - see SARS script
# handle bad filenames separately
class BadPathError(Exception):
    """Exception raised when a path to be created contains illegal characters or reserved
    filenames. Workaround to be able to distinguish from incorrect user input in classic mode.
    """
    pass


# TODO: how long does it take for the run dir to be created?
def find_rundir(run_dir: pathlib.Path, minion_basedir: pathlib.Path) -> pathlib.Path:
    """Check whether run directory exists as full path or directory in MinION dir and
    adjust path of run directory accordingly.
    Arguments:
        run_dir:        absolute or relative path to run directory
        minion_basedir: absolute path to directory of MinION results
    Returns:
        The unchanged run directory if it exists, or the full path to the directory
        within the MinION dir if this was given.
    """
    # we only need to do something if the directory doesn't exist:
    if not run_dir.exists():
        # ...check in MinION dir
        if (minion_basedir / run_dir).exists():
            run_dir = minion_basedir / run_dir
        # if neither of them exists, complain and stop
        else:
            raise FileNotFoundError(f"{str(run_dir)} or {str(minion_basedir / run_dir)} "
                                    f"does not exist \n"
                                    f"Aborting pipeline...")
    # Check if fastq_pass folder exist
    check_fastq_pass = list(run_dir.glob("rawdata/*/fastq_pass"))
    if not check_fastq_pass:
        raise FileNotFoundError(f"fastq_pass folder(s) not found in expected location:\n"
                                f"{str(run_dir)}/rawdata/*/fastq_pass\n"
                                f"Ensure correct directory and/or directory structure is used.\n"
                                f"Aborting pipeline...")

    logger.info(f"Data is retrieved from following folders: \n "
                f"{str([str(fastq_dir) for fastq_dir in check_fastq_pass])}")
    return run_dir


def get_run_name(runsheet: pathlib.Path) -> str:
    """Extract the run name from the runsheet (specified in the column RUNxxxx-INI) # TODO - is it?

    Arguments:
        runsheet:   the path to the runsheet for the run

    Returns:
        The run's name.
    Raises:
        KeyError:   if the runsheet is missing the column for the run name
        ValueError: if the run name hasn't been given in the runsheet
    """
    run_name_col = "RUNxxxx-INI"
    sheet_data = pd.read_excel(runsheet, skiprows = 1, nrows = 2, usecols="A:D")
    if run_name_col not in sheet_data.columns:
        raise KeyError("No column giving the run name found in the runsheet.")
    run_name = sheet_data[run_name_col].squeeze()
    if pd.isna(run_name):
        raise ValueError("No run name given in the runsheet.")
    return run_name


# read runsheet
def process_runsheet(runsheet_path: pathlib.Path,
                     active_config: Dict[str, Any] = workflow_config) -> pd.DataFrame:
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


# get sequencing time from manual input
def ask_seq_time(default_seq_time: float) -> timedelta:
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
        return timedelta(hours=default_seq_time)
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
        return timedelta(hours=sequencing_time)
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
    output_dir_path = get_clean_outdir(output_dir_path)
    logger.info(f"Saving results to {output_dir_path}")
    return output_dir_path


# get the command to run the pipeline
if __name__ == "__main__":
    arg_parser = ArgumentParser(description = "Run the Nanopore amplicon analysis pipeline")
    arg_parser.add_argument("--rundir", help="Full path or name of sequencing folder")
    arg_parser.add_argument("--runsheet", help="Path to runsheet")
    arg_parser.add_argument("--outdir",
                            help="Path to output directory "
                                 "(default: "
                                 f"{workflow_config['paths']['output_base_path']}/[name of rundir])")
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
                                   f" default: {workflow_config['seq_run_duration_hours']})")
    args = arg_parser.parse_args()
    # we assume this is commandline mode unless started without arguments
    manual_mode = False
    if len(sys.argv) == 1:
        manual_mode = True
    # otherwise set up terminal mode
    else:
        # load config if present - we need to do this early since it contains mode information
        if args.workflow_config_file:
            default_config_file = pathlib.Path(args.workflow_config_file).resolve()
            with open(default_config_file, "r", encoding = "utf-8") as config_file:
                workflow_config = yaml.safe_load(config_file)
            # if this is the *only* argument we also change over into manual mode
            # -> three arguments: script name, flag, path
            if len(sys.argv) == 3:
                manual_mode = True
    # load debug settings from config (either the one we loaded or the default)
    debug_run = workflow_config["debug"]
    # load sequencing time settings
    seq_time = timedelta(hours=workflow_config["seq_run_duration_hours"])
    seq_run_fudge_factor = timedelta(hours = 1)
    if manual_mode:
        # lots of typing, so allow tab completion of paths
        readline.set_completer_delims('\t\n=')
        readline.parse_and_bind("tab: complete")
        # ...and greet the user nicely
        print(f"### Nanopore {workflow_config['amplicon_type']} analysis")
        print("# Setup analysis -------------------------------")
        rundir = pathlib.Path(input("Type full path or name of Nanopore "
                                    "sequencing folder and press enter: ").strip().strip("'"))
        seq_time = ask_seq_time(workflow_config["seq_run_duration_hours"])
        runsheet = pathlib.Path(input("Output directory will be based on experiment name."
                                      "\n"
                                      "Enter path to runsheet: ").strip().strip("'")).resolve()
    else:
        # if you're entering this from the commandline you should specify these
        if not args.runsheet:
            raise ValueError("Runsheet not specified.")
        if not args.rundir:
            raise ValueError("Nanopore run directory not specified.")
        runsheet = pathlib.Path(args.runsheet).resolve()
        rundir = pathlib.Path(args.rundir).resolve()
        if args.run_time:
            seq_time = timedelta(hours = args.run_time)

    # check runsheet
    runsheet_data = process_runsheet(runsheet, workflow_config)
    # set up use of LIS features if enabled - TODO: do we only use them for the runsheet check?
    if workflow_config["lab_info_system"]["use_lis_features"]:
        lis_report = workflow_config["lab_info_system"]["lis_report"]
        check_runsheet.check_against_lis(runsheet_data, lis_report, active_config = workflow_config)
    else:
        lis_report = None
    # get output dir from experiment name + amplicon type
    experiment_name = (f"{helpers.extract_nanopore_run_name(runsheet)}"
                       f"-{workflow_config['amplicon_type']}")
    # set output dir
    if args.outdir:  # can only be given in commandline mode
        output_dir = pathlib.Path(args.outdir)
    else:
        output_dir = pathlib.Path(workflow_config["paths"]["output_base_path"]) / experiment_name

    # in manual mode, we'll give the user a chance to fix the path
    if manual_mode:
        ask_output_dir(output_dir)
    # we check the output dir - first check in commandline mode, second in manual mode
    # if something still is broken, or the commandline version has been given a wrong path,
    # we yell at the user and fail
    output_dir = get_clean_outdir(output_dir)
    # we can check for whether this is a test and/or a continued pipeline here - in manual mode the
    # arguments will be false, so it defaults to pipeline defaults
    if args.test_run:
        debug_run = True
        append_to_databases = False
    continue_pipeline = args.continue_pipeline
    # we'll have to handle creating our folders ourselves - catch duplicate dirs here!
    if not output_dir.exists():
        output_dir.mkdir(parents = True)
    else:
        if not continue_pipeline:
            raise FileExistsError("The desired output directory already exists.")
    # create dir for snakemake logs
    if not (output_dir / "logs").exists():
        (output_dir / "logs").mkdir()
    # start pipeline (in Docker container)
    logger.info("Pipeline is now waiting for sequencing to finish...")
    # set up sequencing run
    seq_run = monitor_run.AmpliconRun(sequence_dir = rundir, outdir = output_dir,
                                      runsheet = runsheet,
                                      active_config = workflow_config,
                                      configfile = default_config_file,
                                      test_run = debug_run)
    # calculate waiting time and check interval
    total_time = seq_time + seq_run_fudge_factor
    check_interval = workflow_config["check_interval_seconds"]
    # wait and start - TODO: give pattern more nicely?
    # detach here - keep start log
    if os.fork():
        sys.exit()
    analysis_run = monitor_run.start_on_file_found(seq_run, "*/final_summary*.txt",
                                                   dry_run = args.dry_run,
                                                   watch_timeout = total_time.seconds,
                                                   watch_interval = check_interval)
    logger.info(f"Started pipeline with command {' '.join(analysis_run.args)}")
