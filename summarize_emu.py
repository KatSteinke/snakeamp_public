"""Summarize Emu reports."""

__author__ = "Kat Steinke"

import logging
import math
import pathlib
import re

from argparse import ArgumentParser
from functools import reduce
from typing import Any, Dict, List

import numpy as np
import pandas as pd
import yaml

import helpers
import pipeline_config
import version
__version__ = version.__version__


default_config_file = pipeline_config.default_config_file
workflow_config = pipeline_config.WORKFLOW_DEFAULT_CONF

logger = logging.getLogger("summarize_emu")
logger.setLevel(logging.INFO)
console_log = logging.StreamHandler()
console_log.setLevel(logging.WARNING)
logger.addHandler(console_log)


def report_species_per_barcode(emu_counts: pathlib.Path,
                               active_config: Dict[str, Any] = workflow_config) -> pd.DataFrame:
    """Extract estimated species counts from Emu output (with estimated counts, --keep_counts)
     and recalculate read percentage to include unclassified reads.
     If LIS data is to be used, sample material is added from the LIS report.

    Arguments:
        emu_counts:     the path to Emu's SAMPLE_rel-abundance.tsv file
        active_config:  the config file to use


    Returns:
        Estimated counts and relative abundance for each species, as well as a blank column
        for approval status, for the sample given.

    Raises:
        ValueError: if the name cannot be extracted or if relative abundance does not sum to 1
    """
    # check if name can be extracted to begin with - TODO: nicer flow
    name_and_barcode = None
    name_only = None
    all_names_pattern = helpers.get_id_pattern(active_config['sample_number_settings'][
                                                   'sample_number_format'],
                                               active_config['sample_number_settings'][
                                                   'negative_control'],
                                               active_config['sample_number_settings'][
                                                   'positive_control'])
    sample_name_pattern = re.compile(r"(?P<full_sample_name>"
                                     r"(?P<name_only>"
                                     f"{all_names_pattern.pattern})"
                                     f"_{active_config['barcode_format']})"
                                     r"_rel-abundance\.tsv")
    find_sample_name = re.search(sample_name_pattern, emu_counts.name)
    if find_sample_name:
        sample_name_groups = find_sample_name.groupdict()
        name_and_barcode = sample_name_groups.get("full_sample_name")
        name_only = sample_name_groups.get("name_only")
    if not name_and_barcode:
        raise ValueError(f"File name {emu_counts.name} does not conform to the expected format "
                         "([SAMPLE]_[BARCODE]_rel-abundance.tsv)."
                         " Sample name could not be extracted.")
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
    # reindex so the species stays outside the multiindexed columns
    emu_read_counts = emu_read_counts.set_index("species", drop = True)
    # note down barcode
    barcode_header = [name_and_barcode] * len(emu_read_counts.columns)
    report_headers = [barcode_header]
    # TODO: should we get the translation etc. in a separate function?
    if active_config["lab_info_system"]["use_lis_features"]:
        if active_config["sample_number_settings"]["positive_control"]:  # TODO: move out to separate function
            positive_control_pattern = re.compile("|".join(active_config["sample_number_settings"][
                                                    "positive_control"].keys()))
        else:
            # "unmatchable" regex so nothing gets seen as a positive control when we don't have one
            positive_control_pattern = re.compile('(?!.*)')
        if active_config["sample_number_settings"]["negative_control"]:
            negative_control_pattern = re.compile(active_config["sample_number_settings"][
                                                      "negative_control"])
        else:
            negative_control_pattern = re.compile('(?!.*)')
        sample_material = ""
        # we only want to load the LIS report if we need it
        if not (re.match(positive_control_pattern, name_only)
                or re.match(negative_control_pattern, name_only)):
            prefix_mapping = helpers.get_number_letter_combination(
                active_config["sample_number_settings"][
                    "number_to_letter"],
                active_config["sample_number_settings"][
                    "sample_numbers_in"],
                active_config["sample_number_settings"][
                    "sample_numbers_out"])
            lab_info_data = pd.read_csv(active_config["lab_info_system"]["lis_report"],
                                        encoding = "latin1")
            sample_format_sheet = re.compile(active_config["sample_number_settings"]["format_in_sheet"])
            # start by translating the sample number
            start_pattern = re.compile(r"^" + helpers.parse_out_group_pattern(sample_format_sheet,
                                                                              "sample_type").pattern)
            name_translate = re.sub(start_pattern, lambda match: prefix_mapping.get(match.group(),
                                                                                    match.group()),
                                    name_only)
            component_order_lis = {value: key for key, value in
                                   re.compile(active_config[
                                                  "sample_number_settings"][
                                                  "format_in_lis"]).groupindex.items()}
            name_translate = helpers.rearrange_sample_number(name_translate, sample_format_sheet,
                                                             component_order_lis)
            sample_material = lab_info_data[lab_info_data["prøvenr"] == name_translate]["prøvekategori"].squeeze()
        material_header = [sample_material] * len(emu_read_counts.columns)
        report_headers.append(material_header)

    report_headers.append(emu_read_counts.columns)
    emu_read_counts.columns = pd.MultiIndex.from_arrays(report_headers)
    return emu_read_counts


# combine all on species
def merge_emu(emu_reports: List[pd.DataFrame]) -> pd.DataFrame:
    """Merge Emu reports for all samples. https://stackoverflow.com/a/44338256/15704972

    Arguments:
        emu_reports:   Filtered reports for all samples

    Returns:
        Reports for all samples merged on species ID, sorted by sample name
    """
    combined_report = reduce(lambda left_df, right_df: pd.merge(left_df, right_df,
                                                                how = "outer",
                                                                left_index = True,
                                                                right_index = True),
                             emu_reports)
    combined_report = combined_report.sort_index(level = 0, axis = "columns")
    return combined_report


def merge_all_in_emu_dir(emu_dir: pathlib.Path,
                         active_config: Dict[str, Any] = workflow_config) -> pd.DataFrame:
    """Merge all Emu reports in the supplied directory.

    Arguments:
        emu_dir:        the directory containing all Emu reports
                        (name format: SAMPLE_rel-abundance.tsv)
        active_config:  the config file to use

    Returns:
        All Emu reports in the directory whose names match the name format combined.

    Raises:
        FileNotFoundError:  if the directory does not contain any Emu reports
    """
    emu_reports = list(emu_dir.glob("*_rel-abundance.tsv"))
    if not emu_reports:
        raise FileNotFoundError(f"No Emu reports found in {emu_dir}.")
    all_reports = []
    for emu_report in sorted(emu_reports, key = lambda report: report.name):
        print(emu_report)
        try:
            emu_data = report_species_per_barcode(emu_report, active_config)
        except ValueError as value_err:
            print("Whoops")
            logger.error(f"Error in {emu_report}:\n"
                         f"{value_err}\n"
                         "Empty results will be added to the merged summary.")
            # we know the file matches the pattern
            all_names_pattern = helpers.get_id_pattern(active_config['sample_number_settings'][
                                                           'sample_number_format'],
                                                       active_config['sample_number_settings'][
                                                           'negative_control'],
                                                       active_config['sample_number_settings'][
                                                           'positive_control'])
            sample_name_pattern = re.compile(r"(?P<full_sample_name>"
                                             f"{all_names_pattern.pattern}"
                                             f"_{active_config['barcode_format']})"
                                             r"_rel-abundance\.tsv")
            find_sample_name = re.search(sample_name_pattern,
                                         emu_report.name)
            sample_name_groups = find_sample_name.groupdict()
            sample_name = sample_name_groups.get("full_sample_name")
            emu_data = pd.DataFrame(index = pd.Index(data = ["unassigned"], name = "species"),
                                    columns = pd.MultiIndex.from_arrays([[sample_name,
                                                                          sample_name,
                                                                          sample_name],
                                                                         ["abundance_from_all",
                                                                          "estimated counts",
                                                                          "medtages"]]),
                                    data = [[np.nan, np.nan, ""]])
        all_reports.append(emu_data)

    all_merged = merge_emu(all_reports)
    return all_merged


def write_to_sheets(merged_report: pd.DataFrame, outfile: pathlib.Path) -> None:
    """Export abundance and estimated counts for each sample to one combined and two separate sheets
    in a given output file.

    Arguments:
        merged_report:  the Emu report for all samples
        outfile:        the file to which the reports should be written
    """
    # Pylint complains here but it's a bug
    with (pd.ExcelWriter(path = outfile) as outfile_writer):  # pylint: disable=abstract-class-instantiated
        merged_report.to_excel(outfile_writer, sheet_name = "overview")
        merged_report.loc[:, pd.IndexSlice[:,
                                           ["abundance_from_all"]]].to_excel(outfile_writer,
                                                                             sheet_name = "abundance")
        merged_report.loc[:, pd.IndexSlice[:,
                                           ["estimated counts"]]].to_excel(outfile_writer,
                                                                           sheet_name = "count")


if __name__ == "__main__":
    arg_parser = ArgumentParser(description = "Combine all Emu reports in a given directory")
    arg_parser.add_argument("indir", help = "Directory containing all Emu reports to summarize")
    arg_parser.add_argument("--outfile",
                            help = "File to write Emu results to (default: emu_summarized.xlsx)",
                            default = "emu_summarized.xlsx")
    arg_parser.add_argument("--workflow_config_file",
                            help="Config file for run (overrides default config given in script, "
                                 "can be overridden by commandline options)")
    args = arg_parser.parse_args()
    if args.workflow_config_file:
        default_config_file = pathlib.Path(args.workflow_config_file).resolve()
        with open(default_config_file, "r", encoding = "utf-8") as config_file:
            workflow_config = yaml.safe_load(config_file)
    input_dir = pathlib.Path(args.indir)
    output_file = pathlib.Path(args.outfile)
    merged_emu = merge_all_in_emu_dir(input_dir, active_config = workflow_config)
    write_to_sheets(merged_emu, output_file)
