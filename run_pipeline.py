"""Wrapper script to start the pipeline from the terminal."""

__author__ = "Kat Steinke"

import logging
import pathlib
import re
import subprocess

from argparse import ArgumentParser
from typing import Any, Dict

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
# read runsheet


# check runsheet format
def validate_runsheet_format(run_sheet: pathlib.Path,
                             active_config: Dict[str, Any] = workflow_config) -> None:
    """Identify wrong sample name or barcode formats in runsheet.

    Arguments:
        run_sheet:      Path to runsheet to check
        active_config:  the configuration to use

    Raises:
        ValueError: if sample IDs are malformed or missing
    """
    # check if runsheet was formatted correctly, throw error and print everything that's wrong in one go
    # get only the relevant columns here: KMA nr, barkode; skip the first three lines
    # (not usable for parsing)
    sheet_data = pd.read_excel(run_sheet, usecols="A:C", skiprows=3,
                               dtype={"KMA nr": str, "Barkode NB": str})
    sheet_issues = False
    data_missing = False
    # set up record of issues so they can all be printed at once
    fail_record = "The following issue(s) were detected with the runsheet:"
    # check that sample numbers and barcodes have been entered, fail early if so
    no_sample_ids = sheet_data["KMA nr"].isna().all()
    if no_sample_ids:
        fail_record += "\nNo sample IDs found."
        data_missing = True
    no_barcodes = sheet_data["Barkode NB"].isna().all()
    if no_barcodes:
        fail_record += "\nNo barcodes found."
        data_missing = True
    if data_missing:
        raise ValueError(fail_record)
    # now we can be sure there are sample IDs and barcodes, we can check them
    # start by checking if we have the same amount of sample IDs and barcodes
    amount_sample_ids = sheet_data["KMA nr"].dropna().size
    amount_barcodes = sheet_data["Barkode NB"].dropna().size
    if amount_sample_ids != amount_barcodes:
        sheet_issues = True
        fail_record += f"\nAmount of sample IDs and barcodes don't match. " \
                       f"There are {amount_sample_ids} sample IDs but {amount_barcodes} barcodes."
        # check that positive and negative controls are included
        # for positive controls: see if there are any sample numbers matching the controls
    # check controls if given
    if workflow_config["sample_number_settings"]["positive_control"]:
        positive_control_pattern = "|".join(active_config["sample_number_settings"][
                                                "positive_control"].keys())
        positive_controls_in_sheet = sheet_data["KMA nr"].str.fullmatch(positive_control_pattern,
                                                                        na = False)
        if not positive_controls_in_sheet.any():
            sheet_issues = True
            fail_record += f"No positive controls given in runsheet."
    else:
        positive_control_pattern = ''
    if workflow_config["sample_number_settings"]["negative_control"]:
        negative_control_pattern = active_config["sample_number_settings"]["negative_control"]
        # for negative controls: see if there is anything matching negative control pattern
        negative_controls_in_sheet = sheet_data["KMA nr"].str.fullmatch(negative_control_pattern,
                                                                        na = False)
        if not negative_controls_in_sheet.any():
            sheet_issues = True
            fail_record += "No negative controls given in runsheet."
    else:
        negative_control_pattern = ''

    id_pattern = '^(' \
                 + active_config["sample_number_settings"]["sample_number_format"] \
                 + '|' \
                 + negative_control_pattern \
                 + '|(' \
                 + positive_control_pattern \
                 + '))$'
    fail_ids = sheet_data["KMA nr"][~sheet_data["KMA nr"].apply(str).str.match(id_pattern,
                                                                               na = False)].dropna().tolist()
    if fail_ids:
        sheet_issues = True
        # adapt error message to sample number format - TODO: do we need a second sanity check here
        if active_config['sample_number_settings']['sample_numbers_in'] == "number":
            allowed_start = active_config['sample_number_settings']['number_to_letter'].keys()
        else:
            allowed_start = active_config['sample_number_settings']['number_to_letter'].values()
        # TODO: make sample number length variable?
        fail_record += f"\nSample IDs {fail_ids} are not valid. " \
                       f'Sample IDs must start with {" or ".join(allowed_start)} ' \
                       f'followed by eight numbers (six if leaving out year).' \
                       "Negative controls must be given in the format " \
                       f"{active_config['sample_number_settings']['negative_control']}. " \
                       "Please correct sample IDs in runsheet."


    # check that barcodes have correct format
    fail_barcodes = sheet_data["Barkode NB"][~sheet_data["Barkode NB"].apply(str).str.match(active_config["barcode_format"],
                                                                                            na=False)].dropna().tolist()

    if fail_barcodes:
        sheet_issues = True
        # here we append to the record of issues
        fail_record += f"\nBarcodes {fail_barcodes} are not valid barcodes. " \
                       f"Barcodes must consist of {active_config['barcode_prefix']} " \
                       f"+ a number between 01 and 96."
    # duplicated sample numbers have been checked in the separate runsheet check
    # duplicated barcodes indicate a serious issue though
    if any(sheet_data["Barkode NB"].dropna().duplicated()):
        sheet_issues = True
        duplicated_barcodes = sheet_data["Barkode NB"][sheet_data["Barkode NB"].duplicated()].dropna().unique()
        fail_record += f"\nBarcode(s) {duplicated_barcodes} are duplicated."
    if sheet_issues:
        raise ValueError(fail_record)

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





