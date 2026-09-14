import re
import unittest

from unittest import mock

import pytest

import localization_helpers

class TestGetSupportedLanguages(unittest.TestCase):
    @mock.patch(f'{localization_helpers.__name__}.pathlib.Path.glob')
    def test_fail_no_languages(self, mock_glob):
        """Fail if no languages are found in the locales directory."""
        mock_glob.return_value = []
        error_msg = "No languages found in locales directory."
        with pytest.raises(FileNotFoundError, match=re.escape(error_msg)):
            localization_helpers.get_supported_languages()

    def test_success_find_languages(self):
        """Find the languages for which translations exist."""
        true_supported = {'da', 'en'}
        test_supported = localization_helpers.get_supported_languages()
        assert true_supported == test_supported


class TestSetLanguage(unittest.TestCase):
    def test_fail_unsupported_language(self):
        """Fail if the target language is not supported."""
        error_msg = "tlh is not a supported language. Supported languages are ['da', 'en']."
        with pytest.raises(NotImplementedError, match = re.escape(error_msg)):
            localization_helpers.set_language('tlh')

    def test_success_set_language(self):
        """Successfully set a supported language."""
        _ = localization_helpers.set_language('en')
        expected_translation = "date_received"
        test_translation = _("modtagedato")
        assert expected_translation == test_translation