"""Generic helpers for use in assorted scripts in the pipeline."""

__author__ = "Kat Steinke"

import logging
import pathlib

logger = logging.getLogger("helpers")
logger.setLevel(logging.INFO)
console_log = logging.StreamHandler()
console_log.setLevel(logging.INFO)
logger.addHandler(console_log)


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
