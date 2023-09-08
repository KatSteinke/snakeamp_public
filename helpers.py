"""Generic helpers for use in assorted scripts in the pipeline."""

__author__ = "Kat Steinke"

import logging
import pathlib
import re

from typing import Dict, Union

import pandas as pd
from pandas._libs.missing import NAType

from check_runsheet import logger

logger = logging.getLogger("helpers")
logger.setLevel(logging.INFO)
console_log = logging.StreamHandler()
console_log.setLevel(logging.INFO)
logger.addHandler(console_log)


def get_number_letter_combination(number_to_letter: Dict[str, str], samples_in: str,
                                  samples_out: str) -> Dict[str, str]:
    """Generate correct mapping for a sample number setup from a number to letter mapping.

    Arguments:
        number_to_letter:   The original number to letter mapping
        samples_in:         The format sample numbers have in the runsheet
        samples_out:        The format sample numbers have in the laboratory information system

    Returns:
        A mapping of sample number starts in samples_in format to sample numbers in samples_out
        format.
    """
    # TODO: include some kind of test to ensure this is number to letter?
    if samples_in not in {"number", "letter"}:
        raise ValueError(f"Invalid initial sample format {samples_in}. "
                         f"Sample format can only be number or letter")
    if samples_out not in {"number", "letter"}:
        raise ValueError(f"Invalid desired sample format {samples_out}. "
                         f"Sample format can only be number or letter")
    if samples_in == "number":
        if samples_out == "letter":
            return number_to_letter
        # otherwise, samples_out must be number, so we return a self-to-self mapping
        return {key: key for key in number_to_letter}
    else:  # samples_in is letter
        if samples_out == "number":
            # reverse the dictionary
            return {value: key for (key, value) in number_to_letter.items()}
        # otherwise, samples_out is letter
        return {value: value for value in number_to_letter.values()}


def get_fastq_pass_dir(rundir: pathlib.Path) -> pathlib.Path:
    """Find the fastq_pass directory for the given run directory.

    Arguments:
        rundir: the base directory containing Nanopore sequencing results

    Returns:
        The path to the fastq_pass directory

    Raises:
        FileNotFoundError:  if the fastq_pass directory is not in the expected location
        ValueError:         if there are multiple fastq_pass directories
    """
    # we may need to give the fastq_pass directory directly
    # or a group of dirs in the fastq_pass dir - TODO: do we need to keep this once we move to barcodes?
    if "fastq_pass" in rundir.parts:
        # check if the rundir contains barcodes
        if any((child_dir.name.startswith("barcode") for child_dir in rundir.iterdir())):
            fastq_pass_dir = rundir
        else:
            raise FileNotFoundError("fastq_pass or a subdirectory has been given "
                                    "but no barcode directories were found. "
                                    "Please give the path to the base directory "
                                    "or a directory containing barcode directories ('barcodeXX').")
    else:
        logger.info(f"No barcode directories found in {rundir}.\n"
                    f"Searching for barcodes in {rundir}/rawdata/*/fastq_pass...")
        check_fastq_pass = list(rundir.glob("rawdata/*/fastq_pass"))
        if not check_fastq_pass:
            raise FileNotFoundError(f"fastq_pass folder(s) not found in expected location:\n"
                                    f"{str(rundir)}/rawdata/*/fastq_pass\n"
                                    f"Ensure correct directory and/or directory structure is used.\n"
                                    f"Aborting 16S pipeline...")
        if len(check_fastq_pass) > 1:
            raise ValueError(f"The directory {rundir} contains "
                             f"multiple fastq_pass directories."
                             "Please choose the one containing the fastq files you want to analyze"
                             " and specify the entire path to the fastq_pass directory.")
        fastq_pass_dir = check_fastq_pass[0]
    logger.info(f"Data is retrieved from the following folder:\n"
                f"{fastq_pass_dir}")
    return fastq_pass_dir


def extract_sample_number_part(number_to_check: str, to_extract: str, pattern_in_sheet: re.Pattern,
                               negative_control_pattern: re.Pattern,
                               positive_control_pattern: re.Pattern) -> Union[str, NAType]:
    """Find a part of a sample number described by a match group in the sample format
     used in the runsheet if the sample number is not a control sample.

    Arguments:
        number_to_check:            the sample number from which a component should be extracted
        to_extract:                 the name of the match group describing the component of
                                    the sample number to extract
        pattern_in_sheet:           a regex describing the sample number's elements.
                                    Must contain a match group called sample_type describing
                                    the type prefix format
        negative_control_pattern:   the format in which negative controls are given
        positive_control_pattern:   the format in which positive controls are given

    Returns:
        The component if it could be found; pd.NA otherwise.
    """
    if to_extract not in pattern_in_sheet.groupindex:
        logger.warning(f"Group name {to_extract} not found in sample number pattern."
                       f" Component cannot be extracted.")
        return pd.NA
    # TODO: nicer flow?
    if re.match(negative_control_pattern, number_to_check) \
            or re.match(positive_control_pattern, number_to_check):
        return pd.NA
    if re.match(pattern_in_sheet, number_to_check):
        return re.match(pattern_in_sheet, number_to_check).groupdict().get("sample_type",
                                                                           pd.NA)
    return pd.NA
