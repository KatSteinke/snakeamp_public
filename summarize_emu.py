"""Summarize Emu reports."""

__author__ = "Kat Steinke"

import math
import pathlib
import re

import pandas as pd


# for each Emu report in the directory

def report_species_per_barcode(emu_counts: pathlib.Path) -> pd.DataFrame:
    """Extract estimated species counts from Emu output (with estimated counts, --keep_counts)
     and recalculate read percentage to include unclassified reads.

    Arguments:
        emu_counts: the path to Emu's SAMPLE_rel-abundance.tsv file

    Returns:
        Estimated counts and relative abundance for each species, as well as a blank column
        for approval status, for the sample given.

    Raises:
        ValueError: if the name cannot be extracted or if relative abundance does not sum to 1
    """
    # check if name can be extracted to begin with - TODO: nicer flow
    sample_name = None
    find_sample_name = re.search(r'(?P<sample_name>\w+)_rel-abundance\.tsv', emu_counts.name)
    if find_sample_name:
        sample_name_groups = find_sample_name.groupdict()
        sample_name = sample_name_groups.get("sample_name")
    if not sample_name:
        raise ValueError(f"File name {emu_counts.name} does not conform to the expected format "
                         "([SAMPLE]_rel-abundance.tsv). Sample name could not be extracted.")
    # get read counts per species
    emu_read_counts = pd.read_csv(emu_counts, sep = "\t")
    # check if something is wrong with the abundance as is
    if not math.isclose(sum(emu_read_counts['abundance'].dropna()), 1):
        raise ValueError("Relative abundance does not sum to 1. "
                         "This suggests the result file is broken (missing/extra lines).")
    # recalculate percentage to include unassigned reads
    total_reads = sum(emu_read_counts['estimated counts'])
    emu_read_counts['abundance_from_all'] = emu_read_counts['estimated counts'] / total_reads

    # cut down to required columns and add approval column
    cols_for_report = ["species", "abundance_from_all", "estimated counts", "medtages"]
    emu_read_counts = emu_read_counts.reindex(columns = cols_for_report, fill_value = "")
    # "unassigned" is only noted on the taxid level - fill it in on the species level
    emu_read_counts["species"] = emu_read_counts["species"].fillna(value = "unassigned")
    # note down barcode
    barcode_header = [sample_name] * len(emu_read_counts.columns)
    emu_read_counts.columns = pd.MultiIndex.from_arrays([barcode_header,
                                                         emu_read_counts.columns])
    return emu_read_counts

# combine all on species

