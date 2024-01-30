"""Summarize Emu reports."""

__author__ = "Kat Steinke"

import logging
import math
import pathlib
import re

from argparse import ArgumentParser
from datetime import datetime
from functools import reduce
from typing import Any, Dict, List, NamedTuple

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


def get_lis_information(sample_number: str, lis_report: pd.DataFrame,
                        active_config: Dict[str, Any] = workflow_config) -> pd.DataFrame:
    """Get sample information from LIS (date received, sample category and anatomy).

    Arguments:
        sample_number:  the sample number to look up in the LIS
        lis_report:     a report from the LIS giving sample date ("modtagedato"),
                        category ("prøvekategori") and anatomical location ("anatomi").
        active_config:  the config file to use

    Returns:
        Date received, sample category and anatomical location for the sample (blank for a control).

    Raises:
        KeyError:   if the sample number cannot be found in the LIS report after translation
    """
    (negative_control_pattern,
     positive_control_pattern) = helpers.get_control_patterns(
        active_config["sample_number_settings"]["negative_control"],
        active_config["sample_number_settings"]["positive_control"])
    sample_information = pd.DataFrame(data = {"prøvenr": [sample_number], "modtagedato": [""],
                                              "prøvemateriale": [""],
                                              "anatomi": [""]})
    if not (re.match(positive_control_pattern, sample_number)
            or re.match(negative_control_pattern, sample_number)):
        prefix_mapping = helpers.get_number_letter_combination(
            active_config["sample_number_settings"][
                "number_to_letter"],
            active_config["sample_number_settings"][
                "sample_numbers_output"],
            active_config["sample_number_settings"][
                "sample_numbers_out"])

        sample_format_results = re.compile(active_config["sample_number_settings"]["format_output"])
        sample_format_lis = re.compile(active_config["sample_number_settings"]["format_in_lis"])
        # start by translating the sample number
        name_translate = helpers.translate_sample_number(sample_number, sample_format_results,
                                                         sample_format_lis, prefix_mapping,
                                                         positive_control_pattern,
                                                         negative_control_pattern)
        if name_translate not in lis_report["prøvenr"].tolist():
            error_msg = (f"Sample number {name_translate} (original number: {sample_number}) "
                         "not found in LIS report.")
            raise KeyError(error_msg)
        sample_information = lis_report[lis_report["prøvenr"] == name_translate][["prøvenr",
                                                                                  "modtaget",
                                                                                  "prøvekategori",
                                                                                  "anatomi"]]
        sample_information = sample_information.rename(columns = {"modtaget": "modtagedato",
                                                                  "prøvekategori":
                                                                      "prøvemateriale"})
        sample_information["modtagedato"] = sample_information["modtagedato"].apply(lambda x:
                                                                              datetime.strptime(x,
                                                                                        "%d%m%Y").strftime("%Y-%m-%d"))
        sample_information = sample_information.fillna("").reset_index(drop=True)
    return sample_information


class SampleNameComponents(NamedTuple):
    """Run, sample/isolate number and barcode for a given sample."""
    run_name: str
    sample_name: str
    barcode: str


def extract_name_components(report_name: str, active_config: Dict[str, Any] = workflow_config) \
        -> SampleNameComponents:
    """Extract run name, sample name and barcode from an Emu report's name.

    Arguments:
        report_name:    the name to be parsed
        active_config:  the configuration to be used

    Returns:
        Run name, sample name and barcode encoded in the report's name.
    Raises:
        ValueError: if one or more components are missing
    """
    all_names_pattern = helpers.get_id_pattern(active_config['sample_number_settings'][
                                                   'sample_number_format'],
                                               active_config['sample_number_settings'][
                                                   'negative_control'],
                                               active_config['sample_number_settings'][
                                                   'positive_control'])
    sample_name_pattern = re.compile(r"(?P<run_name>[A-Za-z0-9_æøåÆØÅ-]+)_(?P<name_only>"
                                     f"{all_names_pattern.pattern})"
                                     r"_(?P<barcode>"
                                     f"{active_config['barcode_format']})"
                                     r"_rel-abundance\.tsv")
    find_sample_name = re.search(sample_name_pattern, report_name)
    if find_sample_name:
        sample_name_groups = find_sample_name.groupdict()
        sample_name_components = SampleNameComponents(run_name = sample_name_groups["run_name"],
                                                      sample_name = sample_name_groups["name_only"],
                                                      barcode = sample_name_groups["barcode"])
        return sample_name_components
    raise ValueError(f"File name {report_name} does not conform to the expected format "
                     "([RUN]_[SAMPLE]_[BARCODE]_rel-abundance.tsv)."
                     " Sample name components could not be extracted.")


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
    sample_name_components = extract_name_components(emu_counts.name, active_config)
    # get read counts per species
    emu_read_counts = pd.read_csv(emu_counts, sep = "\t")
    # check if something is wrong with the abundance as is
    if not math.isclose(sum(emu_read_counts['abundance'].dropna()), 1):
        raise ValueError("Relative abundance does not sum to 100%. "
                         "This suggests the result file is broken (missing/extra lines).")
    # recalculate percentage to include unassigned reads
    total_reads = sum(emu_read_counts['estimated counts'])
    # output as percent
    emu_read_counts['abundance_from_all [%]'] = ((emu_read_counts['estimated counts'] / total_reads)
                                                 * 100)
    # round for easier legibility
    emu_read_counts = emu_read_counts.round({"estimated counts": 0, "abundance_from_all [%]": 2})
    emu_read_counts = emu_read_counts.astype({"estimated counts": "Int64"})
    # cut down to required columns and add approval column
    cols_for_report = ["species", "abundance_from_all [%]", "estimated counts", "medtages"]
    emu_read_counts = emu_read_counts.reindex(columns = cols_for_report, fill_value = "")
    # "unassigned" is only noted on the taxid level - fill it in on the species level
    emu_read_counts["species"] = emu_read_counts["species"].fillna(value = "unassigned")
    # deduplicate species names
    # this also sets species as index so we keep it out of the multiindexed columns
    emu_read_counts = emu_read_counts.groupby(by="species").sum()

    # note down relevant information
    run_header = [sample_name_components.run_name] * len(emu_read_counts.columns)
    barcode_header = [sample_name_components.barcode] * len(emu_read_counts.columns)
    name_header = [sample_name_components.sample_name] * len(emu_read_counts.columns)
    report_headers = [run_header, barcode_header, name_header]
    header_names = ["run", "barcode", "prøvenummer"]
    if active_config["lab_info_system"]["use_lis_features"]:
        lis_data = pd.read_csv(active_config["lab_info_system"]["lis_report"],
                               encoding = "latin1", dtype = {"modtaget": str})
        data_from_lis = get_lis_information(sample_name_components.sample_name, lis_data,
                                            active_config)
        # rename sample number if needed - TODO: more prettily!
        name_header = [data_from_lis["prøvenr"].squeeze()] * len(emu_read_counts.columns)
        report_headers[-1] = name_header
        lis_data_cols = ["modtagedato", "prøvemateriale", "anatomi"]
        lis_headers = [[data_from_lis[sample_metadata].squeeze()] * len(emu_read_counts.columns)
                       if pd.notna(data_from_lis[sample_metadata].squeeze())
                       else [""] * len(emu_read_counts.columns)
                       for sample_metadata in lis_data_cols]
        for lis_header in lis_headers:  # TODO: there has to be a prettier solution
            report_headers.append(lis_header)
        header_names.extend(lis_data_cols)
    # add blank PhHV header row
    header_names.append("PhHV")
    report_headers.append([""] * len(emu_read_counts.columns))
    report_headers.append(emu_read_counts.columns)
    header_names += [None]
    emu_read_counts.columns = pd.MultiIndex.from_arrays(report_headers, names = header_names)
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
    # set up fallbacks - sample number is easiest to set up only when we have it..
    fallback_cols = [["", "", ""],
                     ["abundance_from_all [%]",
                      "estimated counts",
                      "medtages"]]
    fallback_names = ["PhHV", None]
    # ... but we don't want to have to check whether we're using LIS features for every sample
    if active_config["lab_info_system"]["use_lis_features"]:
        fallback_cols = [["", "", ""],
                         ["", "", ""],
                         ["", "", ""]] + fallback_cols
        fallback_names = ["modtagedato", "prøvemateriale", "anatomi"] + fallback_names
    for emu_report in sorted(emu_reports, key = lambda reportfile: reportfile.name):
        try:
            emu_data = report_species_per_barcode(emu_report, active_config)
        except ValueError as value_err:
            logger.error(f"Error in {emu_report}:\n"
                         f"{value_err}\n"
                         "Empty results will be added to the merged summary.")
            sample_name_components = extract_name_components(emu_report.name, active_config)
            run_name = sample_name_components.run_name
            sample_name = sample_name_components.sample_name
            barcode = sample_name_components.barcode
            fallback_cols = [[run_name, run_name, run_name],
                             [barcode, barcode, barcode],
                             [sample_name, sample_name, sample_name]] + fallback_cols
            fallback_names = ["run", "barcode", "prøvenummer"] + fallback_names
            fallback_headers = pd.MultiIndex.from_arrays(fallback_cols, names=fallback_names)
            emu_data = pd.DataFrame(index = pd.Index(data = ["unassigned"], name = "species"),
                                    columns = fallback_headers,
                                    data = [[np.nan, np.nan, ""]])

        all_reports.append(emu_data)

    all_merged = merge_emu(all_reports)
    run_names = all_merged.columns.get_level_values("run").unique().tolist()
    if len(run_names) > 1:
        log_msg = f"Data appear to be from multiple runs ({run_names})."
        logger.warning(log_msg)
    return all_merged


def write_to_sheets(merged_report: pd.DataFrame, outfile: pathlib.Path) -> None:
    """Export abundance and estimated counts for each sample to one combined and two separate sheets
    in a given output file.

    Arguments:
        merged_report:  the Emu report for all samples
        outfile:        the file to which the reports should be written
    """
    # Pylint complains here but it's a bug
    with pd.ExcelWriter(path = outfile) as outfile_writer:  # pylint: disable=abstract-class-instantiated
        merged_report.to_excel(outfile_writer, sheet_name = "overview")
        # depending on absence/presence of LIS features we may have more or fewer multiindex levels
        # (sample material/location get added as extra levels)
        # the columns we're interested in are on the last level
        amount_header_cols = merged_report.columns.nlevels - 1
        header_col_slice = [slice(None)] * amount_header_cols
        merged_report.loc[:, (*header_col_slice,
                              "abundance_from_all [%]")].to_excel(outfile_writer,
                                                              sheet_name = "abundance")
        merged_report.loc[:, (*header_col_slice,
                              "estimated counts")].to_excel(outfile_writer,
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
