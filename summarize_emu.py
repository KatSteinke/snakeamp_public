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
    sample_information = pd.DataFrame(data = {"patient": [""],
                                              "prøvenr": [sample_number],
                                              "modtagedato": [""],
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
        sample_information = lis_report[lis_report["prøvenr"] == name_translate][["cprnr.",
                                                                                  "prøvenr",
                                                                                  "modtaget",
                                                                                  "prøvekategori",
                                                                                  "anatomi"]]
        sample_information = sample_information.rename(columns = {"modtaget": "modtagedato",
                                                                  "prøvekategori":
                                                                      "prøvemateriale",
                                                                  "cprnr.": "patient"})
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
        # this might be legit if a fallback file was generated (all reads are unassigned)
        taxids = emu_read_counts['tax_id']
        # TODO: simplify conditions
        if len(taxids) > 1:
            raise ValueError("Relative abundance does not sum to 100%. "
                             "This suggests the result file is broken (missing/extra lines).")
        # if we only have one taxid we can check if it's unassigned
        if not taxids.squeeze() == "unassigned":
            raise ValueError("Relative abundance does not sum to 100%. "
                             "This suggests the result file is broken (missing/extra lines).")
        logger.info(f"All reads for sample {sample_name_components.sample_name} are unassigned.")
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
    emu_read_counts = emu_read_counts.reindex(columns = cols_for_report)
    # "unassigned" is only noted on the taxid level - fill it in on the species level
    emu_read_counts["species"] = emu_read_counts["species"].fillna(value = "unassigned")
    # "medtages" should be a blank string
    emu_read_counts["medtages"] = emu_read_counts["medtages"].fillna(value = "")
    # deduplicate species names
    # this also sets species as index so we keep it out of the multiindexed columns
    emu_read_counts = emu_read_counts.groupby(by="species").sum()
    print(emu_read_counts.to_string())
    # note down relevant information
    run_header = [sample_name_components.run_name] * len(emu_read_counts.columns)
    barcode_header = [sample_name_components.barcode] * len(emu_read_counts.columns)
    name_header = [sample_name_components.sample_name] * len(emu_read_counts.columns)
    report_headers = [run_header, barcode_header, name_header]
    header_names = ["run", "barcode", "prøvenummer"]
    if active_config["lab_info_system"]["use_lis_features"]:
        lis_data = pd.read_csv(active_config["lab_info_system"]["lis_report"],
                               encoding = "latin1", dtype = {"modtaget": str,
                                                             "cprnr.": str})
        data_from_lis = get_lis_information(sample_name_components.sample_name, lis_data,
                                            active_config)
        # rename sample number if needed - TODO: more prettily!
        name_header = [data_from_lis["prøvenr"].squeeze()] * len(emu_read_counts.columns)
        report_headers[-1] = name_header
        lis_data_cols = ["modtagedato", "patient", "prøvemateriale", "anatomi"]
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
    # add blank note row
    header_names.append("notes")
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


def sort_report_samples(emu_report: pd.DataFrame,
                        active_config: Dict[str, Any] = workflow_config) -> pd.DataFrame:
    """Reorder columns in Emu report so that any controls come first and the remaining columns are
    ordered by barcode.

    Arguments:
        emu_report:     multiindexed but unsorted aggregated Emu report
        active_config:  the config file to use

    Returns:
        The Emu report rearranged so that negative and positive controls come first and the rest is
        ordered by barcode (as in the runsheet).
    """
    # get the names of all controls and non-controls and combine them
    (negative_control,
     positive_control) = helpers.get_control_patterns(
        active_config["sample_number_settings"]["negative_control"],
        active_config["sample_number_settings"]["positive_control"])
    controls = re.compile(f"{negative_control.pattern}|{positive_control.pattern}")
    sample_numbers = emu_report.columns.get_level_values("prøvenummer")
    # control columns should be sorted alphabetically to ensure same order
    control_columns = sorted(list({sample_nr
                                   for sample_nr in sample_numbers
                                   if re.search(negative_control, sample_nr)
                                   or re.search(positive_control, sample_nr)}))
    # now get their positions
    control_positions = [position for sample_nr in control_columns
                         for (position, colname) in enumerate(sample_numbers)
                         if colname == sample_nr]
    # non-controls need to be sorted by barcode
    # find all sample numbers not matching control format - predefined slice since it's a lot of writing
    non_control_slice = ~sample_numbers.str.match(controls)
    non_control_header = emu_report.loc[:, non_control_slice].columns.to_frame(index=False)
    # now sort non-control on barcodes since they're guaranteed to be unique
    # (aside from coming in groups of three due to the multiindex setup)
    non_control_barcodes = sorted(list(non_control_header["barcode"].unique()))
    barcodes = emu_report.columns.get_level_values("barcode")
    # get the positions so we can combine them with the controls' positions (based on sample number)
    # TODO: avoid repetition?
    non_control_positions = [position for barcode in non_control_barcodes
                             for (position, colname) in enumerate(barcodes)
                             if colname == barcode]
    reordered_columns_pos = control_positions + non_control_positions
    # reindex with the columns given by positions
    emu_report = emu_report.reindex(pd.MultiIndex.from_tuples([emu_report.columns[column_index]
                                                               for column_index in
                                                               reordered_columns_pos],
                                                              names = emu_report.columns.names),
                                    axis = "columns")
    return emu_report


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
    if not emu_dir.exists():
        raise FileNotFoundError(f"Emu report directory {emu_dir} does not exist")
    if not emu_dir.is_dir():
        raise NotADirectoryError(f"{emu_dir} is a single file. "
                                 f"Please specify the directory containing all Emu reports.")
    emu_reports = list(emu_dir.glob("*_rel-abundance.tsv"))
    if not emu_reports:
        raise FileNotFoundError(f"No Emu reports found in {emu_dir}.")
    all_reports = []
    # set up fallbacks - sample number is easiest to set up only when we have it..
    fallback_cols = [["", "", ""],
                     ["", "", ""],
                     ["abundance_from_all [%]",
                      "estimated counts",
                      "medtages"]]
    fallback_names = ["PhHV", "notes", None]
    # ... but we don't want to have to check whether we're using LIS features for every sample
    if active_config["lab_info_system"]["use_lis_features"]:
        fallback_cols = [["", "", ""],
                         ["", "", ""],
                         ["", "", ""],
                         ["", "", ""]] + fallback_cols
        fallback_names = ["modtagedato", "patient", "prøvemateriale", "anatomi"] + fallback_names
    for emu_report in emu_reports:
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
    # sort here:
    all_merged = sort_report_samples(all_merged, active_config=active_config)
    run_names = all_merged.columns.get_level_values("run").unique().tolist()
    if len(run_names) > 1:
        log_msg = f"Data appear to be from multiple runs ({run_names})."
        logger.warning(log_msg)
    # if we're using LIS data we're giving out a patient ID here
    if "patient" in all_merged.columns.names:
        patient_ids = [pt_id
                       for pt_id in all_merged.columns.get_level_values("patient").unique().tolist()
                       if pt_id]
        # we want to make it clear that this is the ID for this batch of *results*
        # -> use run number from all runs
        # shorter run name - only extract RUNXXXX from names since it's a one-off ID
        # we no longer have zero-padded run numbers
        run_numbers = [re.search(r"run\d{1,4}", run_name, flags = re.IGNORECASE).group(0)
                       if re.search(r"run\d{1,4}", run_name, flags =re.IGNORECASE)
                       else "RUNxxxx"
                       for run_name in run_names]
        all_run_numbers = "_".join(run_numbers)
        patients_to_ids = {patient_id: f"{all_run_numbers}_pt_{patient_index}"
                           for patient_index, patient_id in enumerate(patient_ids)
                           if patient_id}
        all_merged = all_merged.rename(columns = patients_to_ids)
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
        notes = pd.DataFrame(index=pd.Index(merged_report.columns.get_level_values("prøvenummer").unique(),
                                            name="Prøvenummer"),
                             columns = ["notes"])
        notes.to_excel(outfile_writer, sheet_name = "notes")


if __name__ == "__main__":
    arg_parser = ArgumentParser(description = "Combine all Emu reports in a given directory")
    arg_parser.add_argument("indir", help = "Directory containing all Emu reports to summarize")
    arg_parser.add_argument("--outfile",
                            help = "File to write Emu results to (default: emu_summarized.xlsx)",
                            default = "emu_summarized.xlsx")
    arg_parser.add_argument("--outfile_raw",
                            help="File to write raw Emu results to (default: emu_summarized.tsv)",
                            default = "emu_summarized.tsv")
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
    output_file_raw = pathlib.Path(args.outfile_raw)
    merged_emu = merge_all_in_emu_dir(input_dir, active_config = workflow_config)
    write_to_sheets(merged_emu, output_file)
    merged_emu.to_csv(output_file_raw, sep="\t")
