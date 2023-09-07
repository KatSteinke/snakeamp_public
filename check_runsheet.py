"""Check whether the runsheet is correct."""

__author__ = "Kat Steinke"

import logging
import pathlib
import re
import readline
import sys
import warnings

from typing import Any, Dict

import pandas as pd

import helpers
import pipeline_config
import version
__version__ = version.__version__

# import parameters
default_config_file = pipeline_config.default_config_file
workflow_config = pipeline_config.WORKFLOW_DEFAULT_CONF

# convert sample numbers to MADS-friendly format
# see if any weren't found, point at relevant number in runsheet if so

# TODO: refactor to fix code reuse

# TODO: should this raise exceptions or just output an overview?

logger = logging.getLogger("check_runsheet")
logger.setLevel(logging.INFO)
console_log = logging.StreamHandler()
console_log.setLevel(logging.INFO)
plain_messages = logging.Formatter("%(message)s")
console_log.setFormatter(plain_messages)
logger.addHandler(console_log)


def check_by_prefix(sheet_data: pd.DataFrame, lab_data: pd.DataFrame, sheet_prefix: str,
                    lab_data_prefix: str) -> None:
    """Check whether samples of a given sample type (indicated by prefix) are contained in a report
    from the laboratory information system.

    Arguments:
        sheet_data:         Data contained in runsheet
        lab_data:           Data contained in laboratory information system report (minimum:
                            sample numbers and date received)
        sheet_prefix:       Sample number prefix designating desired sample type in sample sheet
        lab_data_prefix:    Sample number prefix corresponding to sheet_prefix in the format used in
                            the LIS report

    Raises:
        ValueError: if samples aren't found in the LIS report
    """
    runsheet_filtered = sheet_data[sheet_data["KMA nr"].str.startswith(sheet_prefix)].copy()
    # check if this leaves us with any data
    if runsheet_filtered.empty:
        raise ValueError(f"No samples with prefix {sheet_prefix} found in runsheet.")
    # piece sample number together: get the prefix, identify the year, then add the last six
    # prefix is already known
    runsheet_filtered["proevenr_prefix"] = sheet_prefix
    # extract sample number
    runsheet_filtered["proevenr_kort"] = runsheet_filtered["KMA nr"].str.slice(start=-6)
    # to find year, check lab information system data and go for date *received*
    # extract relevant samples again
    lab_data_filtered = lab_data[lab_data["prøvenr"].str.startswith(lab_data_prefix)].copy()
    # get bare sample number to match the one from the runsheet
    lab_data_filtered["proevenr_kort"] = lab_data_filtered["prøvenr"].str.slice(start=-6)
    # check if we have duplicates
    if any(lab_data_filtered.duplicated(subset=["proevenr_kort"])):
        raise ValueError("MADS report contains duplicated sample numbers. "
                         "This likely means the report covers multiple years. "
                         "Get a new MADS report with the correct start date.")
    # check if there are mismatches between runsheet and MADS data
    missing_from_mads = runsheet_filtered[~runsheet_filtered["proevenr_kort"].isin(lab_data_filtered["proevenr_kort"])][
        "KMA nr"].dropna().tolist()
    if missing_from_mads:
        raise ValueError(
            f"Samples {missing_from_mads} were not found in MADS report. "
            f"Please check that sample numbers are correct.")


def check_against_lis(sheet_data: pd.DataFrame, lab_report: pathlib.Path,
                      active_config: Dict[str, Any] = workflow_config) -> None:
    """Check if sample numbers are found in laboratory information system report.

    Arguments:
        sheet_data:     DataFrame representation of samples in runsheet
        lab_report:     Path to lab information system's report to check for presence of samples
        active_config:  configuration to use

    Raises:
        ValueError: if sample(s) are not found in LIS report
    """
    # needs prefix parsing, so set this up in the beginning, flip mapping if needed
    prefix_mapping = helpers.get_number_letter_combination(active_config["sample_number_settings"][
                                                               "number_to_letter"],
                                                           active_config["sample_number_settings"][
                                                               "sample_numbers_in"],
                                                           active_config["sample_number_settings"][
                                                               "sample_numbers_out"])
    # load and prepare LIS report

    lab_info_data = pd.read_csv(lab_report, encoding="latin1", dtype = {"afsendt": str,
                                                                        "cprnr.": str,
                                                                        "modtaget": str})
    # sample numbers can be identical except for the prefix
    # -> make subsets of sample sheet and report by prefix
    # remove both negative and positive controls here
    positive_control_pattern = "|".join(active_config["sample_number_settings"][
                                            "positive_control"].keys())
    # extract prefix: numbers or letters
    sample_format_sheet = re.compile(active_config["sample_number_settings"]["format_in_sheet"])
    all_prefixes = sheet_data["KMA nr"].apply(lambda x: re.match(sample_format_sheet,
                                                                 x).groups("sample_type")
                                              if not (re.match(active_config[
                                                                   "sample_number_settings"][
                                                                   "negative_control"],
                                                               x)
                                                      or re.match(positive_control_pattern, x))
                                              else pd.NA)
    unique_prefixes = all_prefixes.dropna().unique()
    # record errors if they happen here
    errors = []
    for prefix in unique_prefixes:
        try:
            check_by_prefix(sheet_data, lab_info_data, prefix, prefix_mapping[prefix])
        except ValueError as value_err:
            errors.append(value_err)

    if errors:
        error_text = "\n".join([str(parsing_error) for parsing_error in errors])
        raise ValueError(f"The following issues were encountered:\n"
                         f"{error_text}")
    else:
        logger.info("The runsheet is correct.")


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
    if active_config["sample_number_settings"]["positive_control"]:
        positive_control_pattern = "|".join(active_config["sample_number_settings"][
                                                "positive_control"].keys())
        positive_controls_in_sheet = sheet_data["KMA nr"].str.fullmatch(positive_control_pattern,
                                                                        na = False)
        if not positive_controls_in_sheet.any():
            sheet_issues = True
            fail_record += f"No positive controls given in runsheet."
    else:
        positive_control_pattern = ''
    if active_config["sample_number_settings"]["negative_control"]:
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
                       f'followed by eight numbers (six if leaving out year). '
        if active_config['sample_number_settings']['negative_control']:
            fail_record += "Negative controls must be given in the format " \
                            f"{negative_control_pattern}. "
        fail_record += "Please correct sample IDs in runsheet."

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


def check_sheet_format(sheet_data: pd.DataFrame, check_barcodes=False,
                       active_config: Dict[str, Any] = workflow_config) -> None:
    """Check if sample numbers are present and in the correct format.

    Arguments:
        sheet_data:     DataFrame representation of samples in runsheet
        check_barcodes: whether or not to check barcodes
        active_config:  configuration to use


    Raises:
        ValueError: if sample IDs are malformed or missing
    """
    sheet_data = sheet_data.copy()
    sheet_issues = False
    data_missing = False
    # set up record of issues so they can all be printed at once
    fail_record = "The following issue(s) were detected with the runsheet:"
    no_sample_ids = sheet_data["KMA nr"].isna().all()
    if no_sample_ids:
        fail_record += "\nNo sample IDs found."
        data_missing = True
    if check_barcodes:
        no_barcodes = sheet_data["Barkode NB"].isna().all()
        if no_barcodes:
            fail_record += "\nNo barcodes found."
            data_missing = True
    if data_missing:
        raise ValueError(fail_record)
    # TODO: just drop NA?
    if check_barcodes:  # TODO - avoid double check?
        # now we can be sure there are sample IDs and barcodes, we can check them
        # start by checking if we have the same amount of sample IDs and barcodes
        amount_sample_ids = sheet_data["KMA nr"].dropna().size
        amount_barcodes = sheet_data["Barkode NB"].dropna().size
        if amount_sample_ids != amount_barcodes:
            sheet_issues = True
            fail_record += f"\nAmount of sample IDs and barcodes don't match. " \
                           f"There are {amount_sample_ids} sample IDs" \
                           f" but {amount_barcodes} barcodes."
            # check that barcodes have correct format
            fail_barcodes = sheet_data["Barkode NB"][
                ~sheet_data["Barkode NB"].apply(str).str.match(active_config["barcode_format"],
                                                               na = False)].dropna().tolist()

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
    # check duplicates early - this only needs to warn, not break
    if any(sheet_data["KMA nr"].dropna().duplicated()):
        duplicated_ids = sheet_data["KMA nr"][sheet_data["KMA nr"].duplicated()].dropna().unique()
        logger.warning(f"Sample number(s) {duplicated_ids} are duplicated. "
                       f"If you are sure you want to sequence the same sample twice, "
                       f"you can ignore this warning.")
    # check controls if given
    if active_config["sample_number_settings"]["positive_control"]:
        positive_control_pattern = "|".join(active_config["sample_number_settings"][
                                                "positive_control"].keys())
        positive_controls_in_sheet = sheet_data["KMA nr"].str.fullmatch(positive_control_pattern,
                                                                        na = False)
        if not positive_controls_in_sheet.any():
            sheet_issues = True
            fail_record += f"\nNo positive controls given in runsheet."
    else:
        positive_control_pattern = ''
    if active_config["sample_number_settings"]["negative_control"]:
        negative_control_pattern = active_config["sample_number_settings"]["negative_control"]
        # for negative controls: see if there is anything matching negative control pattern
        negative_controls_in_sheet = sheet_data["KMA nr"].str.fullmatch(negative_control_pattern,
                                                                        na = False)
        if not negative_controls_in_sheet.any():
            sheet_issues = True
            fail_record += "\nNo negative controls given in runsheet."
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
                                                                               na=False)].dropna().tolist()
    if active_config['sample_number_settings']['sample_numbers_in'] == "number":
        allowed_start = active_config['sample_number_settings']['number_to_letter'].keys()
    else:
        allowed_start = active_config['sample_number_settings']['number_to_letter'].values()
    # TODO: make sample number length variable?
    if fail_ids:
        sheet_issues = True
        fail_record += f"\nSample IDs {fail_ids} are not valid. " \
                       f'Sample IDs must start with {" or ".join(allowed_start)} ' \
                       f'followed by eight numbers (six if leaving out year). '
        if active_config['sample_number_settings']['negative_control']:
            fail_record += "Negative controls must be given in the format " \
                           f"{negative_control_pattern}. "
        fail_record += "Please correct sample IDs in runsheet."
    if sheet_issues:
        raise ValueError(fail_record)


if __name__ == "__main__":
    readline.set_completer_delims('\t\n=')   # allow tab completion of paths
    readline.parse_and_bind("tab: complete")
    # suppress openpyxl warning - not relevant for data processing
    warnings.filterwarnings('ignore',
                            message="Data Validation extension is not supported and will be removed",
                            module="openpyxl")
    print("###Runsheet check")
    run_sheet = pathlib.Path(input("Enter path to runsheet: ").strip().strip("'")).resolve()
    print("Loading runsheet...\n")
    # TODO: handle this part in a function?
    runsheet_data = pd.read_excel(run_sheet, usecols = "A", skiprows = 3,  # don't check CP for now
                                  dtype = {"KMA nr": str})
    runsheet_data = runsheet_data.dropna()
    print("Checking runsheet format....\n")
    # simple error handling, suppressing tracebacks
    try:
        check_sheet_format(runsheet_data)
    except ValueError as value_error:
        logger.error(str(value_error))
        sys.exit()
    lab_info_report = pathlib.Path(workflow_config["lab_info_system"]["lis_report"])
    print("Comparing to samples in MADS......\n")
    try:
        check_against_lis(runsheet_data, lab_info_report)
    except ValueError as value_error:
        logger.error(str(value_error))
        sys.exit()
