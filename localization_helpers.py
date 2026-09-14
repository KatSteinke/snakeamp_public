"""Helper functions for localization"""

__author__ = "Kat Steinke"

import gettext
import pathlib

from typing import Set

def get_supported_languages() -> Set[str]:
    """Identify the languages currently supported by the module.

    Returns:
        The languages currently supported by the module.

    Raises:
        FileNotFoundError:  if no translation files were found
    """
    translation_dirs = (pathlib.Path(__file__).parent / "locales").glob("*/LC_MESSAGES/base.mo")
    # we're looking for the language abbreviation, two steps under the respective "base.mo" file
    supported_languages = {translation_dir.parts[-3] for translation_dir in translation_dirs}
    if not supported_languages:
        raise FileNotFoundError("No languages found in locales directory.")
    return supported_languages


def set_language(target_language: str) -> gettext.NullTranslations.gettext:
    """Set gettext translation for a language if it's part of the set of supported languages.

    Arguments:
        target_language:

    Returns:
        gettext.translation.gettext for the given language if it's supported

    Raises:
         NotImplementedError:   if the target language is not supported
    """
    supported_languages = get_supported_languages()
    target_language = target_language.lower()
    if target_language not in supported_languages:
        raise NotImplementedError(f"{target_language} is not a supported language. "
                                  f"Supported languages are {sorted(list(supported_languages))}.")
    translation = gettext.translation('base',
                                      localedir = pathlib.Path(__file__).parent / "locales",
                                      languages = [target_language])
    translation_gettext = translation.gettext
    return translation_gettext
