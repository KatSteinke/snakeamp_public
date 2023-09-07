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

# TODO: merge with the other runsheet check!
def check_runsheet(sheet_data: pd.DataFrame, lab_report: pathlib.Path, 
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


def check_sample_numbers(sheet_data: pd.DataFrame, 
                         active_config: Dict[str, Any] = workflow_config) -> None:
    """Check if sample numbers are present and in the correct format.

    Arguments:
        sheet_data:     DataFrame representation of samples in runsheet
        active_config:  configuration to use


    Raises:
        ValueError: if sample IDs are malformed or missing
    """
    sheet_data = sheet_data.copy()
    no_sample_ids = sheet_data["KMA nr"].isna().all()
    if no_sample_ids:
        raise ValueError("No sample IDs found.")
    # TODO: just drop NA?
    # check duplicates early - this only needs to warn, not break
    if any(sheet_data["KMA nr"].dropna().duplicated()):
        duplicated_ids = sheet_data["KMA nr"][sheet_data["KMA nr"].duplicated()].dropna().unique()
        logger.warning(f"Sample number(s) {duplicated_ids} are duplicated. "
                       f"If you are sure you want to sequence the same sample twice, "
                       f"you can ignore this warning.")
    # check that positive and negative controls are included
    # for positive controls: see if there are any sample numbers matching the controls
    positive_control_pattern = "|".join(active_config["sample_number_settings"][
                                            "positive_control"].keys())
    positive_controls_in_sheet = sheet_data["KMA nr"].str.fullmatch(positive_control_pattern,
                                                                    na = False)
    if not positive_controls_in_sheet.any():
        raise ValueError(f"No positive controls given in runsheet.")
    # for negative controls: see if there is anything matching negative control pattern
    negative_controls_in_sheet = sheet_data["KMA nr"].str.fullmatch(active_config[
                                                                        "sample_number_settings"][
                                                                        "negative_control"],
                                                                    na = False)
    if not negative_controls_in_sheet.any():
        raise ValueError("No negative controls given in runsheet.")
    id_pattern = '^(' \
                 + active_config["sample_number_settings"]["sample_number_format"] \
                 + '|' \
                 + active_config["sample_number_settings"]["negative_control"] \
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
        raise ValueError(f"\nSample IDs {fail_ids} are not valid. "
                         f"Sample IDs must start with {' or '.join(allowed_start)} "
                         f"followed by eight numbers "
                         f"(six if leaving out year). "
                         f"Negative controls must be given in the format"
                         f" {active_config['sample_number_settings']['negative_control']}. "
                         f"Please correct sample IDs in runsheet.")

# TODO: can't do barcode checks in pre-run runsheet check, but need to in wrapper
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
        check_sample_numbers(runsheet_data)
    except ValueError as value_error:
        logger.error(str(value_error))
        sys.exit()
    lab_info_report = pathlib.Path(workflow_config["lab_info_system"]["lis_report"])
    print("Comparing to samples in MADS......\n")
    try:
        check_runsheet(runsheet_data, lab_info_report)
    except ValueError as value_error:
        logger.error(str(value_error))
        sys.exit()
