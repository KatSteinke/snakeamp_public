"""Wrapper script to start the pipeline from the terminal."""

__author__ = "Kat Steinke"

import logging
import pathlib
import re
import subprocess

from argparse import ArgumentParser

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
# read runsheet

# check runsheet against MADS - TODO: will the year be in the sample number?

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





