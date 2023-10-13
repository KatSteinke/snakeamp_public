"""Check that the pipeline produces the expected output for a test run."""

__author__ = "Kat Steinke"

import logging
import pathlib

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
    expected_organisms = pd.DataFrame(data={"prøvenummer": ["F99123457", "F99123456", "F99123458"],
                                            "organism": ["Staphylococcus agalactiae",
                                                         "Cutibacterium acnes",
                                                         "Lactococcus lactis"]})
    expected_positive_control = pd.DataFrame(data = {"organism": ['Bacillus subtilis',
                                                                  'Staphylococcus aureus',
                                                                  'Listeria monocytogenes',
                                                                  'Salmonella enterica',
                                                                  'Escherichia coli',
                                                                  'Enterococcus faecalis',
                                                                  'Pseudomonas aeruginosa'],
                                                     "abundance": [0.148,
                                                                   0.187,
                                                                   0.037,
                                                                   0.2,
                                                                   0.183,
                                                                   0.124,
                                                                   0.033]})
    results_okay = True
    # get list of tabs - do we have everything
    expected_tabs = {'overview', 'abundance', 'count'}
    with pd.ExcelFile(emu_report) as report_sheet:
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
            expected_headers = pd.DataFrame(data = {"run": ["16S_Run0000-Y20230929-XYZ"] * len(sheet_data.columns),
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
                                                    "prøvenummer": [*["F99123457"] * int(
                                                            len(sheet_data.columns) / num_samples),

                                                                    *["F99123456"] * int(
                                                            len(sheet_data.columns) / num_samples),

                                                                    *["F99123458"] * int(
                                                            len(sheet_data.columns) / num_samples),

                                                                    *["NegK_Sanger"]* int(
                                                            len(sheet_data.columns) / num_samples),

                                                                    *["PosK"] * int(
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
                                                    })
            # are the headers correct? use MultiIndex.to_frame(index=False)
            # strip the "Unnamed" parts out
            sheet_data = sheet_data.rename(columns = lambda colname: "" if "Unnamed" in colname
                                                                       else colname)
            # no need to compare the last column though
            header_cols = sheet_data.columns.to_frame(index = False)[["run",
                                                                      "barcode",
                                                                      "prøvenummer",
                                                                      "modtagedato",
                                                                      "prøvemateriale",
                                                                      "anatomi"]]

            compare_headers = expected_headers.compare(header_cols, result_names = ("expected",
                                                                                    "found"))
            if not compare_headers.empty:
                results_okay = False
                logger.warning("Sample metadata differ from expected sample metadata in tab"
                               f" {sheet}:\n"
                               f"{compare_headers.to_string()}")



    # get the first column for each barcode - this'll be abundance in the overview
    # for the routine samples, is the highest scoring organism what we should expect?
    # for the positive control, are the n highest what we would expect?
    # Is the abundance around where we'd expect it to be?
    return results_okay


