"""Generic helpers for use in assorted scripts in the pipeline."""

__author__ = "Kat Steinke"

#  Copyright (c) 2026 Kat Steinke
#     This program is distributed under version 3 of the GNU General Public License.
#      You should have received a copy of the GNU General Public License
#        along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import logging
import pathlib
import re

from typing import Any, Dict, Optional, Tuple, Union

import pandas as pd
from pandas._libs.missing import NAType

import input_names
import pipeline_config

# import parameters
default_config_file = pipeline_config.default_config_file
workflow_config = pipeline_config.WORKFLOW_DEFAULT_CONF


logger = logging.getLogger("helpers")


sheet_names, lis_names = input_names.load_input_from_config(workflow_config)


class PrettyKeyErrorMessage(str):
    """Workaround to allow formatted explanatory messages when raising a KeyError.
    Taken from https://stackoverflow.com/a/70114007/15704972
    """
    def __repr__(self):
        return str(self)


def get_control_patterns(negative_control: Optional[str] = None,
                         positive_control: Optional[Dict[str, Any]] = None) \
        -> Tuple[re.Pattern, re.Pattern]:
    """Return patterns matching negative and positive controls if given,
     else an unmatchable pattern.

    Arguments:
        negative_control:   the format used for negative controls
        positive_control:   positive controls and their expected results

    Returns:
        Patterns matching negative and positive controls if given, else unmatchable patterns.
    """
    # "unmatchable" regex so nothing gets seen as a control when we don't have one
    negative_control_pattern = re.compile('(?!.*)')
    positive_control_pattern = re.compile('(?!.*)')
    if negative_control:
        negative_control_pattern = re.compile(negative_control)
    if positive_control:
        positive_control_pattern = re.compile("|".join(positive_control.keys()))
    return negative_control_pattern, positive_control_pattern


# TODO: use defaults from config instead?
def get_id_pattern(sample_number_format: str,
                   negative_control: Optional[str]=None,
                   positive_control: Optional[Dict[str, Any]]=None) -> re.Pattern:
    """Get a pattern matching all possible sample numbers and controls, if those are given.

    Argument:
        sample_number_format:   the format used for regular sample numbers
        negative_control:       the format used for negative controls
        positive_control:       positive controls and their expected results

    Returns:
        A pattern that matches any valid sample number or control.
    Raises:
        ValueError: when the sample number pattern is blank.
    """
    if not sample_number_format:
        raise ValueError("Sample number pattern cannot be blank.")
    negative_control_pattern = ""
    if negative_control:
        negative_control_pattern = f"|{negative_control}"
    positive_control_pattern = ""
    if positive_control:
        positive_control_pattern = f"|({'|'.join(positive_control.keys())})"
    pattern_all = re.compile(f"({sample_number_format}"
                             f"{negative_control_pattern}"
                             f"{positive_control_pattern}"
                             f")")
    return pattern_all


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


def check_barcode_dirs(fastq_pass: pathlib.Path) -> None:
    """Check if a fastq_pass directory contains barcode dirs.

    Arguments:
        fastq_pass: the fastq_pass directory to check

    Raises:
        FileNotFoundError:  if there are no barcode dirs
    """
    if not any((child_dir.name.startswith("barcode") for child_dir in fastq_pass.iterdir())):
        raise FileNotFoundError("No barcode directories found in fastq_pass directory.")


def get_fastq_pass_parent(rundir: pathlib.Path) -> pathlib.Path:
    """Find the parent directory of the fastq_pass directory for the given run directory.

    Arguments:
        rundir: the base directory containing Nanopore sequencing results

    Returns:
        The path to the directory containing fastq_pass directory

    Raises:
        FileNotFoundError:  if the fastq_pass directory is not in the expected location
        ValueError:         if there are multiple fastq_pass directories and none has been
                            explicitly specified
    """
    # we may need to give the fastq_pass directory directly
    # or a group of dirs in the fastq_pass dir
    if not rundir.exists():
        raise FileNotFoundError(f"The supplied folder {rundir} does not exist. \n"
                                "Ensure correct directory and/or directory structure is used.\n"
                                "Aborting pipeline...")
    existing_path = rundir.parts
    try:
        # if the fastq_pass directory already is somewhere in the dirs given, use this
        fastq_pass_parts = existing_path[:existing_path.index("fastq_pass") + 1]
        # parts contains the initial "/" - resolve the path to clean this up
        fastq_pass_dir = pathlib.Path("/".join(fastq_pass_parts)).resolve()
        # do we end with something that actually exists?
    except ValueError:
        logger.info(f"Searching for fastq_pass folder in {rundir}...")
        check_fastq_pass = list(rundir.glob("**/fastq_pass"))
        if not check_fastq_pass:
            raise FileNotFoundError(f"fastq_pass folder not found in {rundir} or any subfolders. \n"
                                    "Ensure correct directory and/or directory structure is used.\n"
                                    "Aborting pipeline...")
        if len(check_fastq_pass) > 1:
            raise ValueError(f"The directory {rundir} contains "
                             f"multiple fastq_pass directories."
                             "Please choose the one containing the fastq files you want to analyze"
                             " and specify the entire path to the fastq_pass directory.")
        fastq_pass_dir = check_fastq_pass[0]

    # we're using this at multiple points, at some of which barcode dirs not being present might not
    # be an issue
    try:
        check_barcode_dirs(fastq_pass_dir)
    except FileNotFoundError:
        logger.warning("No barcode directories found in fastq_pass directory. "
                       "This may be due to a delay in copying files from the sequencer, but could"
                       " also mean you have given the wrong path. \n"
                       "Only continue if you are sure. ")
    fastq_pass_parent = fastq_pass_dir.parent
    logger.info(f"Data will be retrieved from the following folder:\n"
                f"{fastq_pass_parent}")
    return fastq_pass_parent


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
        return re.match(pattern_in_sheet, number_to_check).groupdict().get(to_extract,
                                                                           pd.NA)
    return pd.NA


def parse_out_group_pattern(complete_pattern: re.Pattern, group_to_extract: str) -> re.Pattern:
    """Extract the pattern defining a given match group from a pattern containing
     multiple match groups.

    Args:
        complete_pattern:   the pattern from which to extract a match group
        group_to_extract:   the name of the match group to extract

    Returns:
        The pattern defining the match group.
    Raises:
        KeyError:   if the match group is not defined in the pattern
    """
    if group_to_extract not in complete_pattern.groupindex:
        error_msg = f"Group {group_to_extract} not found in named groups."
        raise KeyError(error_msg)
    # regex for extracting: \(\?P\<GROUPNAME\>(?:[^)(]|\((?:[^)(]|\((?:[^)(]|\([^)(]*\))*\))*\))*\)
    # -> get up to three levels of matched parentheses
    named_group_in_pattern = re.compile(r"\(\?P<" \
                                        + group_to_extract \
                                        + r">(?:[^)(]|\((?:[^)(]|\((?:[^)(]|\([^)(]*\))*\))*\))*\)")
    named_group_pattern = re.search(named_group_in_pattern, complete_pattern.pattern)
    if named_group_pattern:
        return re.compile(named_group_pattern.group(0))
    error_msg = f"Group {group_to_extract} not found in named groups."  # TODO: can this even happen now?
    raise KeyError(error_msg)


def rearrange_sample_number(old_sample_number: str, pattern_in: re.Pattern,
                            order_out: Dict[int, str]) -> str:
    """
    Split up a sample number in its components (as given by a regex) and reorder them
    in the order given.

    Arguments:
        old_sample_number:  the sample number in its original order
        pattern_in:         a regex representing the components of the original sample number
        order_out:          the desired order of the components: position to component

    Returns:
        The sample number with components reordered in the desired order.

    Raises:
        KeyError:   if the desired order contains a component not found in the input pattern


    """
    # sanity check if we have everything
    extra_components = set(order_out.values()) - set(pattern_in.groupindex.keys())
    if extra_components:
        error_msg = PrettyKeyErrorMessage(f"Not all desired sample number components could be "
                                          f"found in the original format. "
                                          f"Desired sample number format contains additional "
                                          f"components {extra_components}.")
        raise KeyError(error_msg)
    # find components of the sample number
    sample_components = re.search(pattern_in, old_sample_number)
    # check if any are missing
    if not sample_components:
        raise ValueError(f"No match in sample number {old_sample_number}.")
    # add components in the order given in the input
    sample_reordered = []
    # we're getting this as a key, value tuple - TODO: reverse dict
    component_order = dict(sorted(order_out.items()))
    for component in component_order.values():
        sample_reordered.append(sample_components.group(component))
    # combine components
    new_sample_number = "".join(sample_reordered)
    return new_sample_number


def translate_sample_number(sample_number: str, pattern_in: re.Pattern, pattern_out: re.Pattern,
                            prefix_mapping: Dict[str, str],
                            positive_controls: re.Pattern,
                            negative_controls: re.Pattern) -> str:
    """Translate a sample number by rearranging it to match a desired output pattern and
    optionally substituting prefixes.

    Arguments:
        sample_number:      the original sample number
        pattern_in:         a pattern describing the original format
        pattern_out:        a pattern describing the desired format
        prefix_mapping:     a mapping of what prefix to translate to what
        positive_controls:  format used for positive controls, if any
        negative_controls:  format used for negative controls, if any

    Returns:
        The rearranged and translated sample number.

    Raises:
        KeyError:   if the sample number's prefix is not found in the mapping
        ValueError: if the original sample number does not match the input pattern
                    or the translated sample number does not match the desired output

    """
    # skip controls
    control_patterns = re.compile(f"({positive_controls.pattern})|({negative_controls.pattern})")
    if re.match(control_patterns, sample_number):
        return sample_number
    original_format_match = re.match(pattern_in, sample_number)
    if not original_format_match:
        error_msg = f"Sample number {sample_number} does not match specified input format."
        raise ValueError(error_msg)
    if "sample_type" in pattern_in.groupindex:
        start_pattern = re.compile(r"^" + parse_out_group_pattern(pattern_in,
                                                                  "sample_type").pattern)
        current_start = original_format_match.group("sample_type")
        if current_start not in prefix_mapping:
            error_msg = (f"Prefix {current_start} not found "
                         f"(allowed prefixes are {list(prefix_mapping.keys())}).")
            raise KeyError(error_msg)
        name_translate = re.sub(start_pattern, lambda match: prefix_mapping.get(match.group(),
                                                                                match.group()),
                                sample_number)
    else:
        logger.info("No sample_type given in sample number format specification;"
                    " cannot translate sample type.")
        name_translate = sample_number
    component_order_out = {value: key for key, value in pattern_out.groupindex.items()}
    name_translate = rearrange_sample_number(name_translate, pattern_in,
                                             component_order_out)
    if not re.match(pattern_out, name_translate):
        error_msg = (f"Translated sample number {name_translate} (was {sample_number})"
                     " does not match desired output format.")
        raise ValueError(error_msg)
    return name_translate


def check_experiment_name_problems(experiment_name: str) -> None:
    """Check whether an experiment name contains any parts that may cause issues
     when used as filenames.

    Arguments:
        experiment_name:    the experiment name to check

    Raises:
        ValueError: if the experiment name contains chars that can break something
                    (reserved chars, whitespace, slashes)
    """
    illegal_in_windows = r'[<>:"|?*]'
    if re.search(illegal_in_windows, experiment_name):
        raise ValueError("The run name contains a character that cannot be used in "
                         "Windows filenames.")
    # ...or straight up being a reserved filename
    if pathlib.PureWindowsPath(experiment_name).is_reserved():
        raise ValueError("The run name contains a character that cannot be used in "
                         "Windows filenames.")
    # or it could break something from containing whitespace
    if re.search(r"\s", experiment_name):
        raise ValueError("The run name contains spaces or line breaks.")
    # sometimes people enter multiple initials separated with slashes - this'll create subdirs
    if "/" in experiment_name:
        raise ValueError("The run name contains a forward slash (/). "
                         "This will break the result directory."
                         " Replace forward slashes with underscores (_).")
    logger.debug(f"No issues found with experiment name {experiment_name}.")


def extract_nanopore_run_name(runsheet: pathlib.Path,
                              runsheet_names: input_names.RunsheetNames=sheet_names) -> str:
    """Extract the name of a Nanopore sequencing run from its Excel runsheet.

    Arguments:
         runsheet:          Path to an Excel runsheet containing the Nanopore runsheet
         runsheet_names:    Column names in the runsheet

     Returns:
         The run's name as specified under "RUNxxxx-INI".

    Raises:
        ValueError: if the experiment name contains chars that can break something
                    (reserved chars, whitespace, slashes)
    """
    experiment_sheet = pd.read_excel(runsheet, sheet_name = "Runsheet",
                                   usecols = "A:D", skiprows = 1, nrows=2)
    experiment_name = experiment_sheet.at[0, runsheet_names.experiment_name]
    # the experiment name is used as file names for a lot of things, so catch if it breaks something
    # could break something from containing characters that aren't allowed in Windows
    check_experiment_name_problems(experiment_name)
    return experiment_name

