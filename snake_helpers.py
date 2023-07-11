"""Helper functions for running Snakemake."""

__author__ = "Kat Steinke"

import pathlib

import helpers


def is_gzipped(rundir: pathlib.Path, barcode_number: str) -> bool:
    """Determine whether all files in the barcode directory are gzipped.

    Arguments:
        rundir:         the directory containing sequencing data for the run
        barcode_number: the barcode number for which fastq file format should be identified
                        (as one of either .fastq or .fastq.gz)

    Returns:
        True if fastq files are gzipped, False otherwise
    Raises:
        ValueError: if the extension is unsupported or there are multiple formats in the directory
    """
    fastq_pass = helpers.get_fastq_pass_dir(rundir)
    barcode_dir = fastq_pass / f"barcode{barcode_number}"
    extensions = [read_file.suffix for read_file in barcode_dir.iterdir() if read_file.is_file()]
    if all(extension == ".gz" for extension in extensions):
        return True
    if all(extension == ".fastq" for extension in extensions):
        return False
    # if we haven't returned by now something is wrong
    raise ValueError(f"Extensions {sorted(extensions)} not supported.")

