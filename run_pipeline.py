"""Wrapper script to start the pipeline from the terminal."""

__author__ = "Kat Steinke"

import logging
import pathlib
import re
import subprocess

from argparse import ArgumentParser
from typing import Any, Dict, List, Optional

import pandas as pd

import pipeline_config

# import parameters
default_config_file = pipeline_config.default_config_file
workflow_config = pipeline_config.WORKFLOW_DEFAULT_CONF

# start logging
logger = logging.getLogger("16S_nanopore")
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
                                    f"Aborting ARTIC pipeline...")
    # Check if fastq_pass folder exist
    check_fastq_pass = list(run_dir.glob("rawdata/*/fastq_pass"))
    if not check_fastq_pass:
        raise FileNotFoundError(f"fastq_pass folder(s) not found in expected location:\n"
                                f"{str(run_dir)}/rawdata/*/fastq_pass\n"
                                f"Ensure correct directory and/or directory structure is used.\n"
                                f"Aborting ARTIC pipeline...")

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


# get the command to run the pipeline
def get_pipeline_command(indir: pathlib.Path, outdir: pathlib.Path, runsheet: pathlib.Path,
                         configfile: pathlib.Path = default_config_file,
                         active_config: Dict[str, Any] = workflow_config,
                         debug: Optional[bool] = None) -> List[str]:
    """Generate the command for starting the pipeline.

    Arguments:
        indir:          the directory containing input files for the pipeline
        outdir:         the directory to which results should be output
        runsheet:       the runsheet used for the run
        debug:          whether to run the pipeline in test mode (overrides config setting)
        configfile:     the file containing the configuration for the pipeline
        active_config:  the configuration to use for the pipeline

    Returns:
        The nomad command to start the pipeline
    """
    if debug is None:
        run_as_debug = active_config["debug"]
    else:
        run_as_debug = debug
    nomad_job = "16s-snake-emu-staging" if run_as_debug else "16s-snake-emu-prod"
    nomad_command = ["nomad", "job", "dispatch",
                     "-meta", f"indir={indir}",
                     "-meta", f"outdir={outdir}",
                     "-meta", f"runsheet={runsheet}",
                     nomad_job,
                     configfile]
    return nomad_command


if __name__ == "__main__":
    arg_parser = ArgumentParser(description = "Run the Nanopore 16S analysis pipeline")


