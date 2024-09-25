"""Data structure for setting up input names from config."""

__author__ = "Kat Steinke"

import logging
import pathlib

from typing import Any, Dict, List, NamedTuple, Optional, Tuple, Union

import yaml

import pipeline_config
import version

__version__ = version.__version__

# import parameters
default_config_file = pipeline_config.default_config_file
workflow_config = pipeline_config.WORKFLOW_DEFAULT_CONF

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_log = logging.StreamHandler()
console_log.setLevel(logging.INFO)
logger.addHandler(console_log)


class RunsheetNames(NamedTuple):
    """Column names in the runsheet.

    Attributes:
        sample_number:      column for sample number
        barcode:            column for barcode in Nanopore runsheet
        experiment_name:    column for experiment name
        amplicon_type:      column for amplicon type

    """
    sample_number: str
    barcode: str
    experiment_name: str
    amplicon_type: str


class LISDataNames(NamedTuple):
    """Column names in LIS data.
    Attributes:
      sample_number:    column for sample number in LIS data
      sample_type:      column for sample type (as part of unique sample identifier)
      isolate_number:   column for isolate number (as part of unique isolate identifier)
      patient_id:       column for patient identifier
      date_received:    column for date the sample was received
      material:         column for sample material
      anatomy:          column for anatomical location of the sample
      indication:       column for indication

    """
    sample_number: str
    sample_type: str
    patient_id: str
    date_received: str
    material: str
    anatomy: str
    indication: str
    # isolate number (as part of unique isolate identifier)
    isolate_number: Optional[str] = None


# TODO: should this be configurable (using Nanopore vs Illumina?)
def set_up_runsheet_names(runsheet_config: Dict[str, Union[str, List[str]]]) -> RunsheetNames:
    """Initialize runsheet column names from the runsheet name section of the config.

    Arguments:
        runsheet_config:    the runsheet name section of the config

    Returns:
        Runsheet column and sheet names as a RunsheetNames namedtuple

    Raises:
        KeyError:   if a column name is not filled out
    """
    if not all(runsheet_config.values()):
        missing_names = sorted([column for column, colname in runsheet_config.items()
                                if not colname])
        raise KeyError(f"Column names cannot be blank. Please add name(s) for {missing_names}"
                       " in the run_sheet section of your input settings.")
    runsheet_config = {column: str(colname) for column, colname in runsheet_config.items()}

    current_names = RunsheetNames(**runsheet_config)
    return current_names


# TODO: can we allow for the absence of columns here or is this too much?
def set_up_lis_names(lis_name_config: Dict[str, str]) -> LISDataNames:
    """Initialize LIS column names from the LIS data section of the config.

    Arguments:
        lis_name_config:    the LIS data section of the config

    Returns:
        LIS column names as a LISDataNames namedtuple
    Raises:
        KeyError:   if a column name is not filled out
    """
    optional_cols = {"isolate_number"}
    missing_names = sorted([column for column, colname in lis_name_config.items()
                            if (not colname and column not in optional_cols)])
    if missing_names:
        raise KeyError(f"Column names cannot be blank. Please add name(s) for {missing_names}"
                       " in the lis_columns section of your input settings.")

    # we need to watch out for optional columns here
    lis_name_config = {column: str(lis_name_config[column]) for column in lis_name_config}
    for column in lis_name_config:
        if not lis_name_config[column]:
            lis_name_config[column] = None
    current_names = LISDataNames(**lis_name_config)
    return current_names


def load_input_settings(settings_file: pathlib.Path) \
        -> Tuple[RunsheetNames, LISDataNames]:
    """Load input settings from an input settings YAML file and return settings for each input type.

    Arguments:
        settings_file:   an input settings .yaml file
                         following config/input_settings/input_template.yaml

    Returns:
        RunsheetNames and LISDataNames extracted from the settings file.
    Raises:
         FileNotFoundError:  if the settings file doesn't exist
    """
    if not settings_file.exists():
        raise FileNotFoundError(f"Input config file {settings_file} does not exist.")
    with open(settings_file, "r", encoding = "utf-8") as read_settings:
        all_name_settings = yaml.safe_load(read_settings)
        sheet_names = set_up_runsheet_names(all_name_settings["run_sheet"])
        lis_names = set_up_lis_names(all_name_settings["lis_columns"])
    return sheet_names, lis_names


def load_input_from_config(active_config: Dict[str, Any]) \
        -> Tuple[RunsheetNames, LISDataNames]:
    """Load input settings from the input file given in the workflow config file.
     If the path to the input file is a relative path, it's assumed to be relative to the
      repository root.

    Arguments:
        active_config:  The config to use

    Returns:
        RunsheetNames and LISDataNames extracted from the settings file
         specified in the config.

    """
    input_settings = pathlib.Path(active_config["input_names"])
    if not input_settings.is_absolute():
        input_settings = pathlib.Path(__file__).parent / input_settings
    (sheet_names, lis_names) = load_input_settings(input_settings)
    logger.debug(f"Input settings loaded from {input_settings}")
    return sheet_names, lis_names