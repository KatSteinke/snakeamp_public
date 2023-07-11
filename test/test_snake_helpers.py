import pathlib
import re
import unittest

import pytest

import snake_helpers


class TestIsGzipped(unittest.TestCase):
    rundir = pathlib.Path(__file__).parent / "data" / "snake_helpers" / "test_dir"

    def test_get_fastq(self):
        barcode = "01"
        assert not snake_helpers.is_gzipped(self.rundir, barcode)

    def test_get_gzipped(self):
        barcode = "02"
        assert snake_helpers.is_gzipped(self.rundir, barcode)

    def test_fail_mixed(self):
        barcode = "42"
        error_msg = "Extensions ['.fastq', '.gz'] not supported."
        with pytest.raises(ValueError, match=re.escape(error_msg)):
            snake_helpers.is_gzipped(self.rundir, barcode)