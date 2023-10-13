"""Check that the pipeline produces the expected output for a test run."""

__author__ = "Kat Steinke"

import logging
import pathlib

import pandas as pd

# start logging
logger = logging.getLogger("QATest")
logger.setLevel(logging.INFO)
console_log = logging.StreamHandler()
console_log.setLevel(logging.INFO)
logger.addHandler(console_log)

# check whether emu file has been created to start with
def check_files_present(output_dir: pathlib.Path) -> bool:
    """Check if all expected files are present in the output directory.

    Arguments:
        output_dir: the directory that should contain result files

    Returns:
        True if all files are present, False otherwise.
    """
    # find emu file
    emu_files = list(output_dir.glob("*_emu-combined.xlsx"))
    if not emu_files:
        logger.warning("Emu report is missing. Cannot evaluate Emu results.")
        return False
    if len(emu_files) > 1:
        logger.warning("Multiple Emu summaries found, need only one. "
                       "Cannot evaluate Emu results.")
        return False
    return True



# check whether emu report file contains everything that's needed:
def check_emu_result_file(emu_report: pathlib.Path) -> bool:
    """Check if Emu report contains all data and if the results are correct.

    Arguments:
        emu_report: the path to the Emu report to be checked

    Returns:
        True if all results are correct, False otherwise.
    """
    # get list of tabs - do we have everything
    # load all we have
# are the headers correct? In all the tabs? use MultiIndex.to_frame(index=False)
# for the routine samples, is the highest scoring organism what we should expect?
# for the positive control, are the n highest what we would expect?
# Is the abundance around where we'd expect it to be?
# the same for all the tabs - TODO: turn this into a function?


