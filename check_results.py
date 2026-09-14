"""Check that the pipeline produces the expected output for a test run."""

__author__ = "Kat Steinke"

#  Copyright (c) 2026 Kat Steinke
#     This program is distributed under version 3 of the GNU General Public License.
#      You should have received a copy of the GNU General Public License
#        along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import logging
import math
import pathlib
import sys

from argparse import ArgumentParser
from typing import List

import pandas as pd
import yaml

import pipeline_config
import localization_helpers
import set_log
import version

__version__ = version.__version__

# import parameters
DEFAULT_CONFIG_FILE = pipeline_config.default_config_file
WORKFLOW_CONFIG = pipeline_config.WORKFLOW_DEFAULT_CONF


# start logging
logger = logging.getLogger("QATest")


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
    emu_raw = list(output_dir.glob("*_emu-combined.tsv"))
    if not emu_raw:
        logger.warning("Raw TSV backup of Emu report is missing.")
        return False
    return True


# check whether emu report file contains everything that's needed:
def check_emu_result_file(emu_report: pathlib.Path,
                          report_language: str = WORKFLOW_CONFIG["language"]) -> bool:
    """Check if Emu report contains all data and if the results are correct.

    Arguments:
        emu_report:         the path to the Emu report to be checked
        report_language:    the language of the result file

    Returns:
        True if all results are correct, False otherwise.
    """
    _ = localization_helpers.set_language(target_language = report_language)
    num_samples = 5
    expected_organisms = pd.DataFrame(data={"organism": ["Streptococcus agalactiae",
                                                         "Cutibacterium acnes",
                                                         "Lactococcus lactis"]},
                                      index = pd.Index(["F99123457", "F99123456", "F99123458"],
                                                       name = _("prøvenr")))
    expected_positive_control = pd.DataFrame(data = {"abundance": [15.89,
                                                                   19.58,
                                                                   3.89,
                                                                   19.17,
                                                                   18.14,
                                                                   12.55,
                                                                   5.82,
                                                                   3.03]},
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
                                       header = [0,  # run name
                                                 1,  # version
                                                 2,  # barcode
                                                 3,  # sample number
                                                 4,  # "modtagedato",
                                                 5,  # "patient",
                                                 6,  # "prøvemateriale",
                                                 7,  # "anatomi",
                                                 8,  # "indikation",
                                                 9,  # "total_before_qc",
                                                 10,  # "total_after_qc",
                                                 11,  # "human",
                                                 12,  # PhHV
                                                 13,  # notes
                                                 14])  # abundance/counts/...
            # if it's the overview sheet it'll have a PhHV column, the others don't need one
            # we're not going to compare everything in the PhHV column
            # so don't count this when generating expected data
            amount_compared_cols = len(sheet_data.columns)
            cols_per_sample = int(amount_compared_cols / num_samples)
            # dynamically generate expected headers since some of them might be blank
            expected_headers = pd.DataFrame(data = {"run":
                                                    ["NANO_Amplicon_Y20990101_RUN0001_XYZ-16S"]
                                                        * amount_compared_cols,
                                                    "pipeline_version":
                                                    [f"Version_{__version__}"]
                                                    * amount_compared_cols,
                                                    "barcode": [*["RB62"] * cols_per_sample,
                                                                *["RB64"] * cols_per_sample,
                                                                *["RB31"] * cols_per_sample,
                                                                *["RB51"] * cols_per_sample,
                                                                *["RB60"] * cols_per_sample,
                                                                ],
                                                    _("modtagedato"): [*[""] * cols_per_sample,
                                                                    *[""] * cols_per_sample,
                                                                    *["2021-01-02"] * cols_per_sample,
                                                                    *["2021-01-02"] * cols_per_sample,
                                                                    *["2021-01-02"] * cols_per_sample
                                                                    ],
                                                    _("patient"): [*[""] * cols_per_sample,
                                                                *[""] * cols_per_sample,
                                                                *["RUN0001_pt_0"]
                                                                 * cols_per_sample,
                                                                *["RUN0001_pt_0"]
                                                                 * cols_per_sample,
                                                                *["RUN0001_pt_0"]
                                                                 * cols_per_sample],
                                                    _("prøvemateriale"): [*[""] * cols_per_sample,
                                                        *[""] * cols_per_sample,
                                                        *["Hjerneventrikelvæske <liquor>"]
                                                         * cols_per_sample,
                                                        *["Podning"] * cols_per_sample,
                                                        *["Spinalvæske"] * cols_per_sample,
                                                        ],
                                                    _("anatomi"): [*[""] * cols_per_sample,
                                                                *[""] * cols_per_sample,
                                                                *["Shunt ""(hjerneventrikel)"]
                                                                 * cols_per_sample,
                                                                *["Svælg/tonsil"] * cols_per_sample,
                                                                *[""] * cols_per_sample
                                                                ],
                                                    _("indikation"): [*[""] * cols_per_sample,
                                                                *[""] * cols_per_sample,
                                                                *["!!!"] * cols_per_sample,
                                                                *[""] * cols_per_sample,
                                                                *["en eller anden lang tekst"
                                                                  "<Break/>der ikke kan være på en"
                                                                  " linje i MADS"] * cols_per_sample
                                                                ],
                                                    "total_before_qc": [*[3953] * cols_per_sample,
                                                                         *[167662]
                                                                          * cols_per_sample,
                                                                         *[136886]
                                                                          * cols_per_sample,
                                                                         *[1592] * cols_per_sample,
                                                                        *[21876] * cols_per_sample,
                                                                        ],
                                                    "total_after_qc":
                                                        [*[294] * cols_per_sample,
                                                         *[91561]
                                                          * cols_per_sample,
                                                         *[64426]
                                                          * cols_per_sample,
                                                         *[373] * cols_per_sample,
                                                         *[6831] * cols_per_sample,
                                                         ],
                                                        "human": [*[2] * cols_per_sample,
                                                                  *[6]
                                                                   * cols_per_sample,
                                                                  *[10850]
                                                                   * cols_per_sample,
                                                                  *[0] * cols_per_sample,
                                                                  *[655] * cols_per_sample
                                                                  ]
                                                    }, index = pd.Index([*["NegK_Sanger"]
                                                                          * cols_per_sample,
                                                                    *["PosK"] * cols_per_sample,
                                                                    *["F99123457"]
                                                                     * cols_per_sample,
                                                                    *["F99123456"]
                                                                     * cols_per_sample,
                                                                    *["F99123458"]
                                                                     * cols_per_sample
                                                                    ],
                                                                        name= _("prøvenr")))
            # are the headers correct? use MultiIndex.to_frame(index=False)
            # strip the "Unnamed" parts out
            sheet_data = sheet_data.rename(columns = lambda colname: "" if "Unnamed" in str(colname)
                                                                     else colname)

            # we don't need to compare approval
            header_cols = sheet_data.columns.to_frame(index = False)[["run",
                                                                      "pipeline_version",
                                                                      "barcode",
                                                                      _("prøvenummer"),
                                                                      _("modtagedato"),
                                                                      _("patient"),
                                                                      _("prøvemateriale"),
                                                                      _("anatomi"),
                                                                      _("indikation"),
                                                                      "total_before_qc",
                                                                      "total_after_qc",
                                                                      "human"]]
            # check if we have the same version before we do any more comparison
            report_version = header_cols["pipeline_version"].unique().squeeze()
            if report_version != f"Version_{__version__}":
                raise ValueError(f"Report was created with {report_version}, "
                                 f"is being checked with Version_{__version__}.\n"
                                 "Cannot check reports from a different version of the pipeline.")
            header_cols = header_cols.rename(columns={_("prøvenummer"): "prøvenr"})
            header_cols = header_cols.set_index("prøvenr")
            header_cols[_("modtagedato")] = (pd.to_datetime(header_cols[_("modtagedato")])
                                          .apply(lambda x: x.strftime("%Y-%m-%d")
                                                          if pd.notnull(x)
                                                          else ""))
            try:
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
            except ValueError:
                results_okay = False
                logger.warning("Sample metadata labels differ from expected sample metadata"
                               " - could not compare. \n"
                               f"Expected index: {expected_headers.index}\n"
                               f"Found index: {header_cols.index}\n"
                               f"Expected columns: {expected_headers.columns}\n"
                               f"Found columns: {header_cols.columns}")
            # run extra checks on the overview sheet - TODO: can we avoid checking twice?
            if sheet == "overview":
                # get the first column for each barcode - this'll be abundance in the overview
                # TODO: can we handle the slicing more nicely?
                amount_header_cols = sheet_data.columns.nlevels - 1
                header_col_slice = [slice(None)] * amount_header_cols
                abundances = sheet_data.loc[:, (*header_col_slice, "abundance")]
                # we don't need the extra information now - just keep sample numbers
                abundances.columns = abundances.columns.get_level_values(_("prøvenummer"))
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
                    logger.warning("Incorrect organism for one or more samples. "
                                   "Expected organism(s):\n"
                                   f"{compare_organisms.sort_index().to_string()}")
                # for the positive control, are the n highest what we would expect?
                n_expected_species = len(expected_positive_control.index)
                species_found = (abundances.loc[:,["PosK"]]
                                 .sort_values(by = "PosK", ascending = False)[:n_expected_species])
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
                                                                            abs_tol = 0.5),
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


def check_all_qc(results_dir: pathlib.Path,
                 report_language: str = WORKFLOW_CONFIG["language"]) -> bool:
    """Run all QC checks on a result directory and report success/failure.

    Arguments:
        results_dir:        Directory containing test run results to evaluate
        report_language:    the language of the report to check

    Returns:
        True if all QC checks pass, False otherwise
    """
    files_present = check_files_present(results_dir)
    if not files_present:
        return False
    emu_file = list(results_dir.glob("*_emu-combined.xlsx"))[0]
    check_emu = check_emu_result_file(emu_file, report_language = report_language)
    if not check_emu:
        return False
    logger.info("All QC checks passed")
    return True

def check_results(start_args: List[str]) -> None:
    """Run the result check for the supplied result directory.

    Arguments:
        start_args: the arguments to start the script with

    """
    arg_parser = ArgumentParser(description = "Check whether results of a test run match "
                                              "expected results")
    arg_parser.add_argument("result_dir",
                            help = "Directory containing test run results to evaluate")
    arg_parser.add_argument("-l", "--logfile",
                            help = "File to write log to "
                                   "(default: logs/pipeline_qa.log in result dir)",
                            default = None)
    arg_parser.add_argument("--workflow_config_file", help = "Config file to use",
                            default = None)
    args = arg_parser.parse_args(start_args)
    # TODO make this less messy
    pipeline_logger = set_log.get_stream_log("QATest")
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
    pipeline_logger.addHandler(log_file)
    active_config = WORKFLOW_CONFIG
    if args.workflow_config_file:
        config_file = pathlib.Path(args.workflow_config_file).resolve()
        with open(config_file, "r", encoding = "utf-8") as read_config:
            active_config = yaml.safe_load(read_config)
    check_qc = check_all_qc(result_dir, report_language = active_config["language"])
    if not check_qc:
        logger.error("One or more QC steps failed. Check log for details.")
        set_log.clean_up_handlers(pipeline_logger)
        sys.exit(1)
    set_log.clean_up_handlers(pipeline_logger)
    sys.exit(0)



if __name__ == "__main__":
    check_results(sys.argv[1:])
