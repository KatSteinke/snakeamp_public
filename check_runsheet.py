"""Check whether the runsheet is correct."""

__author__ = "Kat Steinke"

import logging
import pathlib
import re
import readline
import sys
import warnings

from argparse import ArgumentParser
from typing import Any, Dict, List

import pandas as pd
import yaml

import helpers
import pipeline_config
import set_log
import version

from helpers import extract_sample_number_part

__version__ = version.__version__

# import parameters
DEFAULT_CONFIG_FILE = pipeline_config.default_config_file
WORKFLOW_CONFIG = pipeline_config.WORKFLOW_DEFAULT_CONF

# convert sample numbers to MADS-friendly format
# see if any weren't found, point at relevant number in runsheet if so

# TODO: refactor to fix code reuse

# TODO: should this raise exceptions or just output an overview?

logger = logging.getLogger("check_runsheet")
logger.setLevel(logging.DEBUG)


# TODO: sample year is now already spliced in or should be
def check_by_prefix(sheet_data: pd.DataFrame, lab_data: pd.DataFrame, sheet_prefix: str,
                    lab_data_prefix: str,
                    active_config: Dict[str, Any] = WORKFLOW_CONFIG) -> None:
    """Check whether samples of a given sample type (indicated by prefix) are contained in a report
    from the laboratory information system.

    Arguments:
        sheet_data:         Data contained in runsheet
        lab_data:           Data contained in laboratory information system report (minimum:
                            sample numbers and date received)
        sheet_prefix:       Sample number prefix designating desired sample type in sample sheet
        lab_data_prefix:    Sample number prefix corresponding to sheet_prefix in the format used in
                            the LIS report
        active_config:      configuration to use

    Raises:
        ValueError: if samples aren't found in the LIS report
    """
    runsheet_filtered = sheet_data[sheet_data["Prøvenummer"].str.startswith(sheet_prefix)].copy()
    # check if this leaves us with any data
    if runsheet_filtered.empty:
        raise ValueError(f"No samples with prefix {sheet_prefix} found in runsheet.")
    # piece sample number together: get the prefix, identify the year, then add the last six
    # prefix is already known
    runsheet_filtered["proevenr_prefix"] = sheet_prefix
    # extract sample number
    (negative_control_pattern,
     positive_control_pattern) = helpers.get_control_patterns(active_config["sample_number_settings"]["negative_control"],
                                                              active_config["sample_number_settings"]["positive_control"])
    sample_format_sheet = re.compile(active_config["sample_number_settings"]["format_in_sheet"])
    runsheet_filtered["proevenr_kort"] = runsheet_filtered["Prøvenummer"].apply(lambda x:
                                                                           extract_sample_number_part(x,
                                                                                                      "sample_number",
                                                                                                      sample_format_sheet,
                                                                                                      negative_control_pattern,
                                                                                                      positive_control_pattern))
    # to find year, check lab information system data and go for date *received*
    # extract relevant samples again
    lab_data_filtered = lab_data[lab_data["prøvenr"].str.startswith(lab_data_prefix)].copy()
    # get bare sample number to match the one from the runsheet
    sample_format_lis = re.compile(active_config["sample_number_settings"]["format_in_lis"])
    lab_data_filtered["proevenr_kort"] = lab_data_filtered["prøvenr"].apply(lambda x:
                                                                           extract_sample_number_part(x,
                                                                                                      "sample_number",
                                                                                                      sample_format_lis,
                                                                                                      negative_control_pattern,
                                                                                                      positive_control_pattern))
    # check if we have duplicates
    if any(lab_data_filtered.duplicated(subset=["proevenr_kort"])):
        raise ValueError("MADS report contains duplicated sample numbers. "
                         "This likely means the report covers multiple years. "
                         "Get a new MADS report with the correct start date.")
    # check if there are mismatches between runsheet and MADS data
    missing_from_mads = runsheet_filtered[~runsheet_filtered["proevenr_kort"].isin(lab_data_filtered["proevenr_kort"])][
        "Prøvenummer"].dropna().tolist()
    if missing_from_mads:
        raise ValueError(
            f"Samples {missing_from_mads} were not found in MADS report. "
            "Please check that sample numbers are correct.")
    logger.debug(f"All samples with prefix {sheet_prefix} found in LIS.")


def check_against_lis(sheet_data: pd.DataFrame, lab_report: pathlib.Path,
                      active_config: Dict[str, Any] = WORKFLOW_CONFIG) -> None:
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
    # match column names for sample number in sheet and LIS
    sheet_data = sheet_data.rename(columns = {"Prøvenummer": "prøvenr"})
    # remove both negative and positive controls here
    (neg_control_pattern,
     pos_control_pattern) = helpers.get_control_patterns(active_config["sample_number_settings"]["negative_control"],
                                                         active_config["sample_number_settings"][
                                                             "positive_control"])
    # TODO: might be prettier but it's the only way that shuts the warning up
    non_controls = sheet_data[~(sheet_data["prøvenr"].str.fullmatch(pos_control_pattern)
                                | sheet_data["prøvenr"].str.fullmatch(neg_control_pattern))].copy()
    # get the order of components in the sheetvs the  LIS and rearrange accordingly
    sample_format_sheet = re.compile(active_config["sample_number_settings"]["format_in_sheet"])
    sample_format_lis = re.compile(active_config["sample_number_settings"]["format_in_lis"])
    component_order_lis = {value: key for key, value in
                           sample_format_lis.groupindex.items()}
    extra_components = (set(sample_format_sheet.groupindex.keys())
                        - set(component_order_lis.values()))
    if extra_components:
        logger.info(f"Comparing only {list(component_order_lis.values())} to LIS report. "
                    f"Cannot check if {list(extra_components)} component(s) are correct.")
    non_controls["prøvenr_translate"] = non_controls["prøvenr"].apply(lambda x:
                                                                      helpers.translate_sample_number(
                                                                          x,
                                                                          sample_format_sheet,
                                                                          sample_format_lis,
                                                                          prefix_mapping,
                                                                          pos_control_pattern,
                                                                          neg_control_pattern))

    # left join the rest on the LIS report
    samples_in_lis = non_controls.merge(lab_info_data, how = "left",
                                        left_on = "prøvenr_translate", right_on = "prøvenr",
                                        suffixes = ("_sheet", "_lis"))
    # check if anything is missing
    missing_in_mads = samples_in_lis[samples_in_lis["prøvenr_lis"].isna()]["prøvenr_sheet"].tolist()
    if missing_in_mads:
        raise ValueError(
            f"Samples {sorted(missing_in_mads)} were not found in MADS report. "
            "Please check that sample numbers are correct.")
    logger.info("The runsheet is correct.")


def check_sheet_format(sheet_data: pd.DataFrame, check_barcodes=False,
                       active_config: Dict[str, Any] = WORKFLOW_CONFIG) -> None:
    """Check if sample numbers are present and in the correct format.

    Arguments:
        sheet_data:     DataFrame representation of samples in runsheet
        check_barcodes: whether or not to check barcodes
        active_config:  configuration to use


    Raises:
        ValueError: if sample IDs are malformed or missing
    """
    sheet_issues = False
    data_missing = False
    # set up record of issues so they can all be printed at once - TODO: separate data check function?
    fail_record = "The following issue(s) were detected with the runsheet:"
    no_sample_ids = sheet_data["Prøvenummer"].isna().all()
    if no_sample_ids:
        fail_record += "\nNo sample IDs found."
        data_missing = True
    if check_barcodes:
        no_barcodes = sheet_data["Barkode"].isna().all()
        if no_barcodes:
            fail_record += "\nNo barcodes found."
            data_missing = True
    if data_missing:
        raise ValueError(fail_record)
    # TODO: just drop NA?
    if check_barcodes:  # TODO - avoid double check?
        # now we can be sure there are sample IDs and barcodes, we can check them
        # start by checking if we have the same amount of sample IDs and barcodes
        amount_sample_ids = sheet_data["Prøvenummer"].dropna().size
        amount_barcodes = sheet_data["Barkode"].dropna().size
        if amount_sample_ids != amount_barcodes:
            sheet_issues = True
            fail_record += "\nAmount of sample IDs and barcodes don't match. " \
                           f"There are {amount_sample_ids} sample IDs" \
                           f" but {amount_barcodes} barcodes."
        # check that barcodes have correct format
        fail_barcodes = sheet_data["Barkode"][
            ~sheet_data["Barkode"].apply(str).str.match(active_config["barcode_format"],
                                                           na = False)].dropna().tolist()

        if fail_barcodes:
            sheet_issues = True
            # here we append to the record of issues
            fail_record += f"\nBarcodes {fail_barcodes} are not valid barcodes. " \
                           f"Barcodes must consist of {active_config['barcode_prefix']} " \
                           "+ a number between 01 and 96."
        # duplicated sample numbers have been checked in the separate runsheet check
        # duplicated barcodes indicate a serious issue though
        if any(sheet_data["Barkode"].dropna().duplicated()):
            sheet_issues = True
            duplicated_barcodes = sheet_data["Barkode"][sheet_data["Barkode"].duplicated()].dropna().unique()
            fail_record += f"\nBarcode(s) {duplicated_barcodes} are duplicated."
    # check duplicates early - this only needs to warn, not break
    if any(sheet_data["Prøvenummer"].dropna().duplicated()):
        duplicated_ids = sheet_data["Prøvenummer"][sheet_data["Prøvenummer"].duplicated()].dropna().unique()
        logger.warning(f"Sample number(s) {duplicated_ids} are duplicated. "
                       "If you are sure you want to sequence the same sample twice, "
                       "you can ignore this warning.")
    # check controls if given
    if active_config["sample_number_settings"]["positive_control"]:
        positive_control_pattern = "|".join(active_config["sample_number_settings"][
                                                "positive_control"].keys())
        positive_controls_in_sheet = sheet_data["Prøvenummer"].str.fullmatch(positive_control_pattern,
                                                                        na = False)
        if not positive_controls_in_sheet.any():
            sheet_issues = True
            fail_record += f"\nNo positive controls given in runsheet."
    if active_config["sample_number_settings"]["negative_control"]:
        # for negative controls: see if there is anything matching negative control pattern
        negative_controls_in_sheet = sheet_data["Prøvenummer"].str.fullmatch(
            active_config["sample_number_settings"]["negative_control"],
            na = False)
        if not negative_controls_in_sheet.any():
            sheet_issues = True
            fail_record += "\nNo negative controls given in runsheet."

    id_pattern = helpers.get_id_pattern(
        active_config["sample_number_settings"]["sample_number_format"],
        negative_control = active_config["sample_number_settings"]["negative_control"],
        positive_control = active_config["sample_number_settings"]["positive_control"])
    id_pattern = re.compile(f"^{id_pattern.pattern}$")
    fail_ids = sheet_data["Prøvenummer"][~sheet_data["Prøvenummer"].apply(str).str.match(id_pattern,
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
                       'followed by eight numbers (six if leaving out year). '
        if active_config['sample_number_settings']['negative_control']:  # TODO: clean structure
            fail_record += "Negative controls must be given in the format " \
                           f"{active_config['sample_number_settings']['negative_control']}. "
        fail_record += "Please correct sample IDs in runsheet."
    if sheet_issues:
        raise ValueError(fail_record)


def check_runsheet(runsheet: pathlib.Path, check_barcodes: bool = False,
                   active_config: Dict[str, Any] = WORKFLOW_CONFIG) -> bool:
    """Check whether the runsheet format is correct and the samples are present in the LIS report
    if one is used.

    Arguments:
        runsheet:       the path to the runsheet to check
        check_barcodes: whether to check barcodes
        active_config:  the configuration to use

    Returns:
        True if the runsheet format is correct and the samples are present in the LIS report
    Raises:
        ValueError: if the runsheet format is incorrect or samples are missing
    """
    logger.info("Loading runsheet...")
    # TODO: handle this part in a function?
    runsheet_data = pd.read_excel(runsheet, usecols = "A:B", skiprows = 3,  # don't check CP for now
                                  dtype = {"Prøvenummer": str})
    runsheet_data = runsheet_data.dropna(subset = ["Prøvenummer", "Barkode"], how="all")
    logger.info("Checking runsheet format....")
    # simple error handling, suppressing tracebacks
    check_sheet_format(runsheet_data, check_barcodes, active_config=active_config)
    if active_config["lab_info_system"]["use_lis_features"]:
        lab_info_report = pathlib.Path(active_config["lab_info_system"]["lis_report"])
        logger.info("Comparing to samples in MADS......")
        check_against_lis(runsheet_data, lab_info_report, active_config=active_config)
    # if we haven't crashed by now we're fine
    return True


def run_check(input_args: List[Any]) -> None:
    """Run the runsheet check for the supplied runsheet.

    Arguments:
        input_args: the arguments to start the runsheet check with
    """
    parser = ArgumentParser(description = "Check amplicon runsheet")
    parser.add_argument("--runsheet", help = "Path to runsheet to check")
    parser.add_argument("--workflow_config_file",
                        help = f"The config to use (default: {DEFAULT_CONFIG_FILE}).",
                        default = DEFAULT_CONFIG_FILE)
    args = parser.parse_args(input_args)
    # TODO initialize root logger here!
    pipeline_logger = set_log.get_stream_log("check_runsheet", level="DEBUG")
    # suppress openpyxl warning - not relevant for data processing
    warnings.filterwarnings('ignore',
                            message = "Data Validation extension is not supported and will be removed",
                            module = "openpyxl")
    if args.workflow_config_file:
        default_config_file = pathlib.Path(args.workflow_config_file).resolve()
        with open(default_config_file, "r", encoding = "utf-8") as config_file:
            workflow_config = yaml.safe_load(config_file)
    if args.runsheet:
        run_sheet = pathlib.Path(args.runsheet).resolve()
    else:
        # TODO how to test for this nicely
        pipeline_logger.setLevel(logging.INFO)
        plain_messages = logging.Formatter("%(message)s")
        pipeline_logger.handlers[0].setFormatter(plain_messages)
        pipeline_logger.info("###Runsheet check")
        readline.set_completer_delims('\t\n=')  # allow tab completion of paths
        readline.parse_and_bind("tab: complete")

        run_sheet = pathlib.Path(input("Enter path to runsheet: ").strip().strip("'")).resolve()
    try:
        check_runsheet(run_sheet, active_config = workflow_config)
    except ValueError as value_error:
        logger.error(str(value_error))
        set_log.clean_up_handlers(pipeline_logger)
        sys.exit(1)
    set_log.clean_up_handlers(pipeline_logger)
    sys.exit(0)


if __name__ == "__main__":
    run_check(sys.argv[1:])
