"""Wrapper script to start the pipeline from the terminal."""

__author__ = "Kat Steinke"

import logging
import pathlib
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




