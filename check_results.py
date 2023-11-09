"""Check that the pipeline produces the expected output for a test run."""

__author__ = "Kat Steinke"

import logging
import math
import pathlib
import sys

from argparse import ArgumentParser

import pandas as pd

# start logging
logger = logging.getLogger("QATest")
logger.setLevel(logging.INFO)
console_log = logging.StreamHandler()
console_log.setLevel(logging.INFO)
logger.addHandler(console_log)


# check whether emu file has been created to start with
def check_files_present(output_dir: pathlib.Path) -> bool:
    """Check if all expected files are present in the output directory.

    Arguments:
        output_dir: the directory that should contain result files

    Returns:
        True if all files are present, False otherwise.
    """
    # find emu file
    emu_files = list(output_dir.glob("*_emu-combined.xlsx"))
    if not emu_files:
        logger.warning("Emu report is missing. Cannot evaluate Emu results.")
        return False
    if len(emu_files) > 1:
        logger.warning("Multiple Emu summaries found, need only one. "
                       "Cannot evaluate Emu results.")
        return False
    return True


# check whether emu report file contains everything that's needed:
def check_emu_result_file(emu_report: pathlib.Path) -> bool:
    """Check if Emu report contains all data and if the results are correct.

    Arguments:
        emu_report: the path to the Emu report to be checked

    Returns:
        True if all results are correct, False otherwise.
    """
    num_samples = 5
    expected_organisms = pd.DataFrame(data={"organism": ["Streptococcus agalactiae",
                                                         "unassigned",
                                                         "Lactococcus lactis"]},
                                      index = pd.Index(["F99123457", "F99123456", "F99123458"],
                                                       name = "prøvenr"))
    expected_positive_control = pd.DataFrame(data = {"abundance": [14.80,
                                                                   18.70,
                                                                   3.70,
                                                                   20.00,
                                                                   18.30,
                                                                   12.40,
                                                                   5.80,
                                                                   3.30]},
                                             index = pd.Index(['Bacillus subtilis',
                                                               'Staphylococcus aureus',
                                                               'Listeria monocytogenes',
                                                               'Salmonella enterica',
                                                               'Escherichia coli',
                                                               'Enterococcus faecalis',
                                                               'Limosilactobacillus fermentum',
                                                               'Pseudomonas aeruginosa'],
                                                              name = "organism"))
    results_okay = True
    # get list of tabs - do we have everything
    expected_tabs = {'overview', 'abundance', 'count'}
    with (pd.ExcelFile(emu_report) as report_sheet):
        sheets_in_report = set(report_sheet.sheet_names)
        tabs_found = expected_tabs.intersection(sheets_in_report)
        if len(tabs_found) < len(expected_tabs):
            results_okay = False
            missing_tabs = expected_tabs - tabs_found
            logger.warning(f"Tab(s) {sorted(list(missing_tabs))} not found in Emu report. "
                           "Cannot evaluate results for these tabs.")
        # for each sheet:
        for sheet in tabs_found:
            sheet_data = pd.read_excel(report_sheet, sheet_name = sheet, index_col = 0,
                                       header = [0, 1, 2, 3, 4, 5, 6])
            # dynamically generate expected headers since some of them might be blank
            expected_headers = pd.DataFrame(data = {"run":
                                                        ["16S_Run0000-Y20230929-XYZ"] * len(sheet_data.columns),
                                                    "barcode": [*["RB31"] * int(
                                                            len(sheet_data.columns) / num_samples),
                                                                *["RB51"] * int(
                                                                        len(sheet_data.columns) / num_samples),
                                                                *["RB60"] * int(
                                                                        len(sheet_data.columns) / num_samples),

                                                                *["RB62"] * int(
                                                                        len(sheet_data.columns) / num_samples),

                                                                *["RB64"] * int(
                                                                        len(sheet_data.columns) / num_samples)],
                                                    "modtagedato": [*["2021-01-02"] * int(
                                                            len(sheet_data.columns) / num_samples),

                                                                    *["2021-01-02"]* int(
                                                            len(sheet_data.columns) / num_samples),

                                                                    *["2021-01-02"]* int(
                                                            len(sheet_data.columns) / num_samples),

                                                                   *[""]* int(
                                                            len(sheet_data.columns) / num_samples),

                                                                    *[""]* int(
                                                            len(sheet_data.columns) / num_samples)],
                                                    "prøvemateriale": [
                                                        *["Hjerneventrikelvæske <liquor>"]* int(
                                                            len(sheet_data.columns) / num_samples),

                                                        *["Podning"]* int(
                                                            len(sheet_data.columns) / num_samples),

                                                        *["Spinalvæske"]* int(
                                                            len(sheet_data.columns) / num_samples),

                                                        *[""]* int(
                                                            len(sheet_data.columns) / num_samples),

                                                        *[""]* int(
                                                            len(sheet_data.columns) / num_samples)],
                                                    "anatomi": [*["Shunt (hjerneventrikel)"]* int(
                                                            len(sheet_data.columns) / num_samples),

                                                                *["Svælg/tonsil"]* int(
                                                            len(sheet_data.columns) / num_samples),

                                                                *[""]* int(
                                                            len(sheet_data.columns) / num_samples),

                                                                *[""]* int(
                                                            len(sheet_data.columns) / num_samples),

                                                                *[""]* int(
                                                            len(sheet_data.columns) / num_samples)]
                                                    }, index = pd.Index([*["F99123457"] * int(
                                                            len(sheet_data.columns) / num_samples),

                                                                    *["F99123456"] * int(
                                                            len(sheet_data.columns) / num_samples),

                                                                    *["F99123458"] * int(
                                                            len(sheet_data.columns) / num_samples),

                                                                    *["NegK_Sanger"]* int(
                                                            len(sheet_data.columns) / num_samples),

                                                                    *["PosK"] * int(
                                                            len(sheet_data.columns) / num_samples)],
                                                                        name="prøvenr"))
            # are the headers correct? use MultiIndex.to_frame(index=False)
            # strip the "Unnamed" parts out
            sheet_data = sheet_data.rename(columns = lambda colname: "" if "Unnamed" in str(colname)
                                                                     else colname)
            # no need to compare the last column though
            header_cols = sheet_data.columns.to_frame(index = False)[["run",
                                                                      "barcode",
                                                                      "prøvenummer",
                                                                      "modtagedato",
                                                                      "prøvemateriale",
                                                                      "anatomi"]]
            header_cols = header_cols.rename(columns={"prøvenummer": "prøvenr"})
            header_cols = header_cols.set_index("prøvenr")
            header_cols["modtagedato"] = pd.to_datetime(header_cols["modtagedato"]).apply(lambda x:
                                                                                          x.strftime(
                                                                                              "%Y-%m-%d")
                                                                                          if pd.notnull(x)
                                                                                          else "")
            compare_headers = expected_headers.compare(header_cols, result_names = ("expected",
                                                                                    "found"))
            if not compare_headers.empty:
                results_okay = False
                # for a "proper" header in the overview tab we'll have duplicated entries
                # but we can't assume that so we only deduplicate now
                compare_headers = compare_headers.drop_duplicates()
                logger.warning("Sample metadata differ from expected sample metadata in tab"
                               f" {sheet}:\n"
                               f"{compare_headers.to_string()}")
        if "overview" in sheets_in_report:  # TODO: handle more nicely - avoid having to reload
            overview_sheet = pd.read_excel(report_sheet, sheet_name = "overview", index_col = 0,
                                           header = [0, 1, 2, 3, 4, 5, 6])
            # get the first column for each barcode - this'll be abundance in the overview
            # TODO: can we handle the slicing more nicely?
            amount_header_cols = overview_sheet.columns.nlevels - 1
            header_col_slice = [slice(None)] * amount_header_cols
            abundances = overview_sheet.loc[:, (*header_col_slice, "abundance_from_all")]
            # we don't need the extra information now - just keep sample numbers
            abundances.columns = abundances.columns.get_level_values("prøvenummer")
            # for the routine samples, is the highest scoring organism what we should expect?
            routine_orgs = abundances.loc[:, ["F99123457",
                                              "F99123456",
                                              "F99123458"]]
            # get the index (=name) of the organism with the highest value
            found_organisms = pd.DataFrame(data=routine_orgs.idxmax().rename("organism"))
            found_organisms.index.names = ["prøvenr"]
            compare_organisms = expected_organisms.compare(found_organisms,
                                                           result_names = ("expected",
                                                                           "found"))
            if not compare_organisms.empty:
                results_okay = False
                logger.warning("Incorrect organism for one or more samples. Expected organism(s):\n"
                               f"{compare_organisms.sort_index().to_string()}")
            # for the positive control, are the n highest what we would expect?
            n_expected_species = len(expected_positive_control.index)
            species_found = abundances.loc[:,
                                           ["PosK"]].sort_values(by = "PosK",
                                                                 ascending = False)[:n_expected_species]
            species_found.index.names = ["organism"]
            species_found.columns.name = None
            species_found = species_found.rename(columns={"PosK": "abundance"})
            missing_species = set(expected_positive_control.index) - set(species_found.index)
            extra_species = set(species_found.index) - set(expected_positive_control.index)
            if missing_species or extra_species:
                results_okay = False
                logger.warning("Positive control should contain"
                               f" {sorted(expected_positive_control.index.tolist())}, "
                               f"contains {sorted(species_found.index.tolist())} "
                               f"(missing: {missing_species}, extra: {extra_species}")
            # Is the abundance around where we'd expect it to be?
            abundances_to_compare = expected_positive_control.merge(species_found, how="inner",
                                                                    left_index = True,
                                                                    right_index = True,
                                                                    suffixes=("_expected",
                                                                              "_found"))
            abundances_match = abundances_to_compare.apply(lambda df:
                                                           math.isclose(df["abundance_expected"],
                                                                        df["abundance_found"],
                                                                        rel_tol = 0.001,
                                                                        abs_tol = 0.1),
                                                           axis=1)
            if not all(abundances_match):
                results_okay = False
                abundance_diff = abundances_to_compare[~abundances_match]
                abundance_diff = abundance_diff.rename(columns = lambda colname:
                                                                 str.replace(colname,
                                                                             "abundance_",
                                                                             ""))
                logger.warning("Different abundance in positive control for "
                               f"{abundance_diff.index.tolist()}."
                               " Expected abundance:\n"
                               f"{abundance_diff.to_string()}")
    return results_okay

if __name__ == "__main__":
    arg_parser = ArgumentParser(description = "Check whether results of a test run match "
                                              "expected results")
    arg_parser.add_argument("result_dir", help="Directory containing test run results to evaluate")
    arg_parser.add_argument("-l", "--logfile", help="File to write log to "
                                                    "(default: logs/pipeline_qa.log in result dir)",
                            default = None)
    args = arg_parser.parse_args()
    result_dir = pathlib.Path(args.result_dir)
    if args.logfile:
        logfile_path = pathlib.Path(args.logfile)
    else:
        logfile_path = result_dir / "logs" / "pipeline_qa.log"
    # log to file
    log_file = logging.FileHandler(logfile_path)
    log_file.setLevel(logging.INFO)
    logfile_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    log_file.setFormatter(logfile_formatter)
    logger.addHandler(log_file)

    # TODO: better way to log pass/fail?
    qc_passes = []
    files_present = check_files_present(result_dir)
    qc_passes.append(files_present)
    if files_present:
        emu_file = list(result_dir.glob("*_emu-combined.xlsx"))[0]
        check_emu = check_emu_result_file(emu_file)
        qc_passes.append(check_emu)
    if not all(qc_passes):
        logger.error("One or more QC steps failed. Check log for details.")
        sys.exit(1)
    sys.exit(0)

