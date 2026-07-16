import logging
import pathlib
import re
import unittest

import pytest

import input_names
from input_names import RunsheetNames, LISDataNames


class TestGetRunsheetNames(unittest.TestCase):
    def test_fail_missing_names(self):
        """Complain if there are blank column or sheet names."""
        error_msg = ("Column names cannot be blank. Please add name(s) for ['barcode',"
                     " 'sample_number'] in the run_sheet section of your input settings.")
        fail_config = {"sample_number": "",
                       "barcode": "",
                       "experiment_name": "Experiment name",
                       "amplicon_type": "Analysis"
                       }
        with pytest.raises(KeyError, match=re.escape(error_msg)):
            input_names.set_up_runsheet_names(fail_config)

    def test_success(self):
        """Successfully read the settings"""
        test_config = {"sample_number": "Prøvenr",
                       "barcode": "Barkode",
                       "experiment_name": "Experiment name",
                       "amplicon_type": "Analysis"
                       }
        expected_result = input_names.RunsheetNames(sample_number = "Prøvenr",
                                                    barcode = "Barkode",
                                                    experiment_name = "Experiment name",
                                                    amplicon_type = "Analysis")
        test_result = input_names.set_up_runsheet_names(test_config)
        assert test_result == expected_result

    def test_success_handle_numbers(self):
        """Handle numbers instead of names for columns - force everything to strings."""
        test_config = {"sample_number": 1,
                       "barcode": "Barkode",
                       "experiment_name": "Experiment name",
                       "amplicon_type": "Analysis"
                       }
        expected_result = input_names.RunsheetNames(sample_number = "1",
                                                    barcode = "Barkode",
                                                    experiment_name = "Experiment name",
                                                    amplicon_type = "Analysis")
        test_result = input_names.set_up_runsheet_names(test_config)
        assert test_result == expected_result


class TestSetupLISData(unittest.TestCase):
    def test_fail_missing_names(self):
        """Fail if column names are missing."""
        fail_msg = ("Column names cannot be blank. Please add name(s) for ['sample_number',"
                    " 'sample_type'] in the lis_columns section of your input settings.")
        fail_config = {"sample_number": "",
                       # sample type (as part of unique sample identifier)
                       "sample_type": "",
                       # isolate number (as part of unique isolate identifier)
                       "isolate_number": "BAKT_NR",
                       "patient_id": "cprnr",
                       "date_received": "modtaget",
                       "material": "prøvekategori",
                       "anatomy": "anatomi",
                       "indication": "Indikation"}
        with pytest.raises(KeyError, match=re.escape(fail_msg)):
            input_names.set_up_lis_names(fail_config)

    def test_success(self):
        """Successfully read the config."""
        test_config = {"sample_number": "PRV_NR",
                       # sample type (as part of unique sample identifier)
                       "sample_type": "P_TYPE",
                       # isolate number (as part of unique isolate identifier)
                       "isolate_number": "BAKT_NR",
                       "patient_id": "cprnr",
                       "date_received": "modtaget",
                       "material": "prøvekategori",
                       "anatomy": "anatomi",
                       "indication": "Indikation"}
        expected_result = input_names.LISDataNames(sample_number = "PRV_NR", sample_type = "P_TYPE",
                                                   patient_id = "cprnr",
                                                   isolate_number = "BAKT_NR",
                                                   date_received = "modtaget",
                                                   material = "prøvekategori",
                                                   anatomy = "anatomi",
                                                   indication = "Indikation")
        test_result = input_names.set_up_lis_names(test_config)
        assert test_result == expected_result

    def test_success_no_isolate(self):
        """Successfully read the config, handling the absence of an optional column."""
        test_config = {"sample_number": "PRV_NR",
                       # sample type (as part of unique sample identifier)
                       "sample_type": "P_TYPE",
                       # isolate number (as part of unique isolate identifier)
                       "isolate_number": "",
                       "patient_id": "cprnr",
                       "date_received": "modtaget",
                       "material": "prøvekategori",
                       "anatomy": "anatomi",
                       "indication": "Indikation"}
        expected_result = input_names.LISDataNames(sample_number = "PRV_NR", sample_type = "P_TYPE",
                                                   patient_id = "cprnr",
                                                   isolate_number = None,
                                                   date_received = "modtaget",
                                                   material = "prøvekategori",
                                                   anatomy = "anatomi",
                                                   indication = "Indikation")
        test_result = input_names.set_up_lis_names(test_config)
        assert test_result == expected_result


    def test_success_numbers(self):
        """Handle columns having numeric names."""
        test_config = {"sample_number": 1,
                       # sample type (as part of unique sample identifier)
                       "sample_type": "P_TYPE",
                       # isolate number (as part of unique isolate identifier)
                       "isolate_number": "BAKT_NR",
                       "patient_id": "cprnr",
                       "date_received": "modtaget",
                       "material": "prøvekategori",
                       "anatomy": "anatomi",
                       "indication": "Indikation"}
        expected_result = input_names.LISDataNames(sample_number = "1", sample_type = "P_TYPE",
                                                   patient_id = "cprnr",
                                                   isolate_number = "BAKT_NR",
                                                   date_received = "modtaget",
                                                   material = "prøvekategori",
                                                   anatomy = "anatomi",
                                                   indication = "Indikation")
        test_result = input_names.set_up_lis_names(test_config)
        assert test_result == expected_result


class TestLoadInputConfig(unittest.TestCase):
    def test_fail_missing_file(self):
        """Fail if the input config file is missing."""
        fail_file = pathlib.Path(__file__).parent / 'data' / 'input_names' / 'no_such_config.yaml'
        error_msg = f"Input config file {fail_file} does not exist."
        with pytest.raises(FileNotFoundError, match=re.escape(error_msg)):
            input_names.load_input_settings(fail_file)

    def test_success_load_input_config(self):
        """Successfully load the different config sections."""
        input_file = pathlib.Path(__file__).parent / 'data' / 'input_names' / 'input_en.yaml'
        sheet_names_en = RunsheetNames(sample_number = "Sample_number",
                                       barcode = "Barcode", experiment_name = "Experiment name",
                                       amplicon_type = "Analysis")
        lis_names_en = LISDataNames(sample_number = "SAMPLENR", sample_type = "SAMPLETYPE",
                                    patient_id = "patient_id",
                                    isolate_number = "BACT_NR",
                                    date_received = "received",
                                    material = "sample_category",
                                    anatomy = "anatomy",
                                    indication = "indication")

        (sheet_names_test,
         lis_names_test) = input_names.load_input_settings(input_file)
        assert sheet_names_en == sheet_names_test
        assert lis_names_en == lis_names_test


class TestLoadFromConfig(unittest.TestCase):
    @pytest.fixture(autouse = True)
    def inject_fixtures(self, caplog):
        self._caplog = caplog

    def test_success_load_from_abs_path(self):
        """Load input settings from an absolute path."""
        test_config = {"input_names": pathlib.Path(__file__).parent / "data"
                                         / "input_names" / "input_en.yaml"}
        log_msg = f"Input settings loaded from {test_config['input_names']}"
        sheet_names_en = RunsheetNames(sample_number = "Sample_number",
                                       barcode = "Barcode",
                                       experiment_name = "Experiment name",
                                       amplicon_type = "Analysis")
        lis_names_en = LISDataNames(sample_number = "SAMPLENR", sample_type = "SAMPLETYPE",
                                    patient_id = "patient_id",
                                    isolate_number = "BACT_NR",
                                    date_received = "received",
                                    material = "sample_category",
                                    anatomy = "anatomy",
                                    indication = "indication")
        with self._caplog.at_level(logger="input_names", level=logging.DEBUG):
            (sheet_names_test,
             lis_names_test) = input_names.load_input_from_config(test_config)
            assert ("input_names", logging.DEBUG, log_msg) in self._caplog.record_tuples
        assert sheet_names_en == sheet_names_test
        assert lis_names_en == lis_names_test

    def test_success_load_from_relative_path(self):
        """Load input settings from a relative path."""
        test_config = {"input_names": "config/input_settings/input_en.yaml"}
        expected_path = (pathlib.Path(__file__).parent.parent / "config" / "input_settings"
                         / "input_en.yaml")
        log_msg = f"Input settings loaded from {expected_path}"
        sheet_names_en = RunsheetNames(sample_number = "Sample_number",
                                       barcode = "Barcode",
                                       experiment_name = "Experiment name",
                                       amplicon_type = "Analysis")
        lis_names_en = LISDataNames(sample_number = "SAMPLENR", sample_type = "SAMPLETYPE",
                                    patient_id = "patient_id",
                                    isolate_number = "BACT_NR",
                                    date_received = "received",
                                    material = "sample_category",
                                    anatomy = "anatomical_location",
                                    indication = "indication")
        with self._caplog.at_level(logger = "input_names", level = logging.DEBUG):
            (sheet_names_test,
             lis_names_test) = input_names.load_input_from_config(test_config)
            assert ("input_names", logging.DEBUG, log_msg) in self._caplog.record_tuples
        assert sheet_names_en == sheet_names_test
        assert lis_names_en == lis_names_test


